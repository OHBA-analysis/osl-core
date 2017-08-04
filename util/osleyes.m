classdef osleyes < handle
	
	% TODO - Add support for different volumes
	% TODO - Add in some text controls

	properties
		clims = {}; % Values below range are transparent!
		colormaps= {}; % Can be name of a colormap function on the path
		current_point = [1 1 1 1]; % Includes time index
	end

	properties(SetAccess=protected)
		niifiles
		colormap_resolution = 255; % Number of colors in each colormap (if not specified as matrix)
	end

	properties(GetAccess=private,SetAccess=private)
		fig
		ax
		h_crosshair
		nii
		coord
		h_img = {} % Cell array of img handles associated with each nii file (there are 3)
		colormap_matrices = {}; % Cached colormaps
		under_construction = true; % Don't render anything while this is true
		motion_active = 0; % Index of which axis to compare against
		lims = nan(3,2); % Axis limits for each MNI dimension (used in plots)
	end

	properties(Dependent)
		valid % Object is valid only if its bound figure is still open
	end

	methods

		function self = osleyes(niifiles)

			if nargin < 1 || isempty(niifiles) 
				niifiles = {fullfile(osldir,'std_masks/MNI152_T1_8mm_brain.nii.gz')};
            end
			
            if ~iscell(niifiles)
                niifiles = {niifiles};
            end

            
			self.fig = figure('Units','normalized');
			set(self.fig,'CloseRequestFcn',@(a,b) delete(self));

			self.ax(1) = axes('Position',[0 0.2 0.3 0.8]);
			self.ax(2) = axes('Position',[0.35 0.2 0.3 0.8]);
			self.ax(3) = axes('Position',[0.7 0.2 0.3 0.8]);
			hold(self.ax(1),'on');
			hold(self.ax(2),'on');
			hold(self.ax(3),'on');
			set(self.fig,'Units','pixels'); % So dragging works by comparing against getpixelposition(ax(j))

			self.nii = []
			self.niifiles = niifiles;
			for j = 1:length(self.niifiles)
				self.nii{j} = load_untouch_nii(self.niifiles{j}); % Do not apply xform/qform
				self.nii{j}.img = double(self.nii{j}.img);
				self.coord{j} = get_coords(self.nii{j}.hdr);
				self.h_img{j}(1) = image(self.coord{j}.y,self.coord{j}.z,permute(self.nii{j}.img(1,:,:,1),[2 3 1])','Parent',self.ax(1),'HitTest','off');
				self.h_img{j}(2) = image(self.coord{j}.x,self.coord{j}.z,permute(self.nii{j}.img(:,1,:,1),[1 3 2])','Parent',self.ax(2),'HitTest','off');
				self.h_img{j}(3) = image(self.coord{j}.x,self.coord{j}.y,permute(self.nii{j}.img(:,:,1,1),[1 2 3])','Parent',self.ax(3),'HitTest','off');
				self.clims{j} = [min(self.nii{j}.img(:)),max(self.nii{j}.img(:))];
				self.colormaps{j} = 'bone';
			end

			self.lims(1,:) = [min(cellfun(@(x) min(x.x),self.coord)) max(cellfun(@(x) max(x.x),self.coord))];
			self.lims(2,:) = [min(cellfun(@(x) min(x.y),self.coord)) max(cellfun(@(x) max(x.y),self.coord))];
			self.lims(3,:) = [min(cellfun(@(x) min(x.z),self.coord)) max(cellfun(@(x) max(x.z),self.coord))];
			set(self.ax(1),'XLim',self.lims(2,:),'YLim',self.lims(3,:));
			set(self.ax(2),'XLim',self.lims(1,:),'YLim',self.lims(3,:));
			set(self.ax(3),'XLim',self.lims(1,:),'YLim',self.lims(2,:));

			set(self.ax(1),'Color','k','View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal');
			set(self.ax(2),'Color','k','View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
			set(self.ax(3),'Color','k','View', [0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');

			self.h_crosshair(1) = plot(self.ax(1),NaN,NaN,'g','HitTest','off');
			self.h_crosshair(2) = plot(self.ax(2),NaN,NaN,'g','HitTest','off');
			self.h_crosshair(3) = plot(self.ax(3),NaN,NaN,'g','HitTest','off');
			set(self.ax,'ButtonDownFcn',@(a,b) axis_click(a,b,self))

			set(self.fig,'WindowButtonDownFcn',@(~,~) self.activate_motion())
			set(self.fig,'WindowButtonMotionFcn',@(a,b) move_mouse(self))
			set(self.fig,'WindowButtonUpFcn',@(~,~) self.deactivate_motion())

			self.under_construction = false; % Enable rendering
			self.refresh_colors();

		end

		function set.colormaps(self,val)
			if self.under_construction
				self.colormaps = val;
				return;
			end

			assert(iscell(val),'Colormaps should be a cell array')
			assert(length(val) == length(self.nii),'Number of colormaps must match number of images');
			for j = 1:length(val)
				if iscell(val{j})
					assert(length(val{j})==2,'Colormap must either be a value or a cell array of length 2');
				end
			end
			self.colormaps = val;
			self.refresh_colors;
		end

		function set.clims(self,val)
			if self.under_construction
				self.clims = val;
				return;
			end

			assert(iscell(val),'clims should be a cell array')
			assert(length(val) == length(self.nii),'Number of clims must match number of images');
			for j = 1:length(self.colormaps)
				assert(length(val{j})==2,'Colour limits must have two elements')
				assert(val{j}(2)>=val{j}(1),'Colour limits must be in ascending order')
				if iscell(self.colormaps{j})
					assert(val{j}(1)>=0,'If using bidirectional limits, colour range must start >= 0')
				end
			end
			self.clims = val;
			self.refresh_colors;
		end

		function set.current_point(self,val)
			% Validate limits
			xl = get(self.ax(3),'XLim');
			yl = get(self.ax(1),'XLim');
			zl = get(self.ax(1),'YLim');

			if val(1) > xl(2) || val(1) < xl(1)
				return
			end

			if val(2) > yl(2) || val(2) < yl(1)
				return
			end

			if val(3) > zl(2) || val(3) < zl(1)
				return
			end

			if length(val) == 3
				val(4) = 1;
			end

			self.current_point = val;
			refresh_slices(self)
		end

		function v = get.valid(self)
			v = ishandle(self.fig);
		end

		function delete(self)
			delete(self.fig);
		end

	end

	methods(Access=private)
		function refresh_slices(self)
			% Given crosshair position, update the slices

			p = self.current_point; % Current point in 3D

			set(self.h_crosshair(1),'XData',[self.lims(2,:) NaN p(2) p(2)],'YData',[p(3) p(3) NaN self.lims(3,:)]);
			set(self.h_crosshair(2),'XData',[self.lims(1,:) NaN p(1) p(1)],'YData',[p(3) p(3) NaN self.lims(3,:)]);
			set(self.h_crosshair(3),'XData',[self.lims(1,:) NaN p(1) p(1)],'YData',[p(2) p(2) NaN self.lims(2,:)]);

			% Now update each slice

			for j = 1:length(self.nii)
				[~,idx(1)] = min(abs(self.coord{j}.x-p(1)));
				[~,idx(2)] = min(abs(self.coord{j}.y-p(2)));
				[~,idx(3)] = min(abs(self.coord{j}.z-p(3)));

				% These are the slice maps - need to convert them to color values now
				d1 = permute(self.nii{j}.img(idx(1),:,:,p(4)),[2 3 1])';
				d2 = permute(self.nii{j}.img(:,idx(2),:,p(4)),[1 3 2])';
				d3 = permute(self.nii{j}.img(:,:,idx(3),p(4)),[1 2 3])';



				set(self.h_img{j}(1),'CData',map_colors(d1,self.colormap_matrices{j},self.clims{j}),'AlphaData',+(abs(d1)>self.clims{j}(1)));
				set(self.h_img{j}(2),'CData',map_colors(d2,self.colormap_matrices{j},self.clims{j}),'AlphaData',+(abs(d2)>self.clims{j}(1)));
				set(self.h_img{j}(3),'CData',map_colors(d3,self.colormap_matrices{j},self.clims{j}),'AlphaData',+(abs(d3)>self.clims{j}(1)));
			end
		end

		function refresh_colors(self)
			% Turn the colormap strings into colormap matrices
			for j = 1:length(self.colormaps)
				if iscell(self.colormaps{j})
                    self.colormap_matrices{j} = cell(2,1);
					self.colormap_matrices{j}{1} = feval(self.colormaps{j}{1},self.colormap_resolution);
					self.colormap_matrices{j}{2} = feval(self.colormaps{j}{2},self.colormap_resolution);
				else
					self.colormap_matrices{j} = feval(self.colormaps{j},self.colormap_resolution);
				end
			end
			self.refresh_slices;
		end

		function activate_motion(self)
			if is_within(get(self.fig,'CurrentPoint'),getpixelposition(self.ax(1)))
				self.motion_active = 1;
			elseif is_within(get(self.fig,'CurrentPoint'),getpixelposition(self.ax(2)))
				self.motion_active = 2;
			elseif is_within(get(self.fig,'CurrentPoint'),getpixelposition(self.ax(3)))
				self.motion_active = 3;
			end
		end

		function deactivate_motion(self);
			self.motion_active = 0;
		end

	end

end

function rgb = map_colors(x,cmap,clim)
	% x - data
	% cmap - color matrix OR cell with two color matrices
	% clim - color limit for this plot

	if iscell(cmap) % We need to treat the positive and negative parts separately
		% These are the positive indices
		pos_idx = min(size(cmap{1},1),round((size(cmap{1},1)-1)*(x-clim(1))/(clim(2)-clim(1)))+1);
        pos_idx(pos_idx<1) = 1;
		neg_idx = min(size(cmap{2},1),round((size(cmap{2},1)-1)*(-x-clim(1))/(clim(2)-clim(1)))+1);
        neg_idx(neg_idx<1) = 1;

		r = zeros(size(x));
		g = zeros(size(x));
		b = zeros(size(x));
		r(x>=0) = cmap{1}(pos_idx(x>=0),1);
		g(x>=0) = cmap{1}(pos_idx(x>=0),2);
		b(x>=0) = cmap{1}(pos_idx(x>=0),3);
		r(x<0) = cmap{2}(neg_idx(x<0),1);
		g(x<0) = cmap{2}(neg_idx(x<0),2);
		b(x<0) = cmap{2}(neg_idx(x<0),3);
	else
		idx = min(size(cmap,1),round((size(cmap,1)-1)*(x-clim(1))/(clim(2)-clim(1)))+1);
		idx(idx<1) = 1;

		r = reshape(cmap(idx,1),size(x));
		g = reshape(cmap(idx,2),size(x));
		b = reshape(cmap(idx,3),size(x));
	end

	rgb = cat(3,r,g,b);

end

function axis_click(a,b,self)
	p = get(a,'CurrentPoint');
	if a == self.ax(1)
		self.current_point(2:3) = p(1,1:2);
	elseif a == self.ax(2)
		self.current_point([1,3]) = p(1,1:2);
	else
		self.current_point(1:2) = p(1,1:2);
	end
end

function c = get_coords(hdr)
	% See https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html#ref0

	i = 0:hdr.dime.dim(2)-1;
	j = 0:hdr.dime.dim(3)-1;
	k = 0:hdr.dime.dim(4)-1;

	c = struct;
	if hdr.hist.sform_code > 0 || (hdr.hist.sform_code == 0 && hdr.hist.qform_code == 0)
		% Use SFORM
		c.x = hdr.hist.srow_x(1) * i + hdr.hist.srow_x(4);
		c.y = hdr.hist.srow_y(2) * j + hdr.hist.srow_y(4);
		c.z = hdr.hist.srow_z(3) * k + hdr.hist.srow_z(4);
	else
		warning('Using qform, not fully tested!')
		% Use QFORM
		a = 0;
		b = hdr.hist.quatern_b;
		c = hdr.hist.quatern_c;
		d = hdr.hist.quatern_d;

		R = [a^2+b^2-c^2-d^2  , 2*b*c-2*a*d     , 2*b*d+2*a*c;
		     2*b*c+2*a*d      , a^2-b^2+c^2-d^2 , 2*c*d-2*a*b;
		     2*b*d-2*a*c      , 2*c*d+2*a*b     , a^2-b^2-c^2+d^2;];

		c.x = R*hdr.dime.pixdim(2)*bsxfun(@times,i(:),[1 0 0])' + hdr.hist.qoffset_x;
		c.y = R*hdr.dime.pixdim(3)*bsxfun(@times,j(:),[0 1 0])' + hdr.hist.qoffset_y;
		c.z = R*hdr.dime.pixdim(1)*hdr.dime.pixdim(4)*bsxfun(@times,k(:),[0 0 1])' + hdr.hist.qoffset_z;
	end

end

function move_mouse(self)
	if isMultipleCall();  return;  end
	if ~self.motion_active; return; end

	if self.motion_active == 1 && is_within(get(self.fig,'CurrentPoint'),getpixelposition(self.ax(1)))
		p = get(self.ax(1),'CurrentPoint');
		self.current_point(2:3) = p(1,1:2);
	elseif self.motion_active == 2 && is_within(get(self.fig,'CurrentPoint'),getpixelposition(self.ax(2)))
		p = get(self.ax(2),'CurrentPoint');
		self.current_point([1,3]) = p(1,1:2);
	elseif self.motion_active == 3 && is_within(get(self.fig,'CurrentPoint'),getpixelposition(self.ax(3)))
		p = get(self.ax(3),'CurrentPoint');
		self.current_point(1:2) = p(1,1:2);
	end
end

function within = is_within(C,pos)
	if C(1) < pos(1) || C(1) > pos(1)+pos(3) || C(2) < pos(2) || C(2) > pos(2)+pos(4)
		within = false;
	else
		within = true;
	end
end

function flag=isMultipleCall()
	% Based on Yair Altman's function
	flag = false; 
	% Get the stack
	s = dbstack();
	if numel(s)<=2
		% Stack too short for a multiple call
		return
	end

	% How many calls to the calling function are in the stack?
	names = {s(:).name};
	TF = strcmp(s(2).name,names);
	count = sum(TF);
	if count>1
		% More than 1
		flag = true; 
	end
end

