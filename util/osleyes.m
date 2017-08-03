classdef osleyes < handle

	properties
		clims = {}; % Values below range are transparent!
		colormaps= {}; % Can be name of a colormap function on the path
		current_point = [1 1 1 1]; % Includes time index
	end

	properties(SetAccess=protected)
		fig
		ax
		nii
		coord
		h_img = {} % Cell array of img handles associated with each nii file (there are 3)
		h_crosshair
		colormap_resolution = 255; % Number of colors in each colormap (if not specified as matrix)
	end

	properties(GetAccess=private,SetAccess=private)
		colormap_matrices = {}; % Cached colormaps 
	end


	methods

		function self = osleyes(niifiles)

			niifiles = {'std_masks/MNI152_T1_8mm_brain.nii.gz'}

			self.fig = figure;
			self.ax(1) = axes();
			self.ax(2) = axes();
			self.ax(3) = axes();
			hold(self.ax(1),'on');
			hold(self.ax(2),'on');
			hold(self.ax(3),'on');
			set([self.fig self.ax],'Units','normalized');

			self.nii = []
			for j = 1:length(niifiles)
				self.nii{j} = load_untouch_nii(niifiles{j}); % Do not apply xform/qform
				self.nii{j}.img = double(self.nii{j}.img);
				self.coord{j} = get_coords(self.nii{j}.hdr);
				self.h_img{j}(1) = image(self.coord{j}.y,self.coord{j}.z,permute(self.nii{j}.img(1,:,:,1),[2 3 1])','Parent',self.ax(1),'HitTest','off');
				self.h_img{j}(2) = image(self.coord{j}.x,self.coord{j}.z,permute(self.nii{j}.img(:,1,:,1),[1 3 2])','Parent',self.ax(2),'HitTest','off');
				self.h_img{j}(3) = image(self.coord{j}.x,self.coord{j}.y,permute(self.nii{j}.img(:,:,1,1),[1 2 3])','Parent',self.ax(3),'HitTest','off');
				self.clims{j} = [min(self.nii{j}.img(:)),max(self.nii{j}.img(:))];
				self.colormaps{j} = 'bone';
			end

			xlims = [min(cellfun(@(x) min(x.x),self.coord)) max(cellfun(@(x) max(x.x),self.coord))];
			ylims = [min(cellfun(@(x) min(x.y),self.coord)) max(cellfun(@(x) max(x.y),self.coord))];
			zlims = [min(cellfun(@(x) min(x.z),self.coord)) max(cellfun(@(x) max(x.z),self.coord))];
			set(self.ax(1),'XLim',ylims,'YLim',zlims);
			set(self.ax(2),'XLim',xlims,'YLim',zlims);
			set(self.ax(3),'XLim',xlims,'YLim',ylims);

			set(self.ax(1),'Color','k','Position',[0 0.2 0.3 0.8],'View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal');
			set(self.ax(2),'Color','k','Position',[0.35 0.2 0.3 0.8],'View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
			set(self.ax(3),'Color','k','Position',[0.7 0.2 0.3 0.8],'View', [0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');

			self.h_crosshair(1) = plot(self.ax(1),NaN,NaN,'g','HitTest','off');
			self.h_crosshair(2) = plot(self.ax(2),NaN,NaN,'g','HitTest','off');
			self.h_crosshair(3) = plot(self.ax(3),NaN,NaN,'g','HitTest','off');
			set(self.ax,'ButtonDownFcn',@(a,b) axis_click(a,b,self))


			self.refresh_colors();

		end

		function set.colormaps(self,val)
			assert(length(val) == length(self.clims),'Number of colormaps must match number of clims');
			for j = 1:length(val)
				if iscell(val{j})
					assert(length(val{j})==2,'Colormap must either be a value or a cell array of length 2');
				end
			end
			self.colormaps = val;
			self.refresh_colors;
		end

		function set.clims(self,val)
			assert(length(val) == length(self.colormaps),'Number of clims must match number of colormaps');
			for j = 1:length(self.colormaps)
				assert(length(val{j}==2),'Colour limits must have two elements')
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

			val(1) = max(xl(1), min(xl(2), val(1)));
			val(2) = max(yl(1), min(yl(2), val(2)));
			val(3) = max(zl(1), min(zl(2), val(3)));

			self.current_point = val;
			refresh_slices(self)
		end


	end

	methods(Access=private)
		function refresh_slices(self)
			% Given crosshair position, update the slices
			p = self.current_point; % Current point in 3D

			set(self.h_crosshair(1),'XData',[get(self.ax(1),'XLim') NaN p(2) p(2)],'YData',[p(3) p(3) NaN get(self.ax(1),'YLim')]);
			set(self.h_crosshair(2),'XData',[get(self.ax(2),'XLim') NaN p(1) p(1)],'YData',[p(3) p(3) NaN get(self.ax(2),'YLim')]);
			set(self.h_crosshair(3),'XData',[get(self.ax(3),'XLim') NaN p(1) p(1)],'YData',[p(2) p(2) NaN get(self.ax(3),'YLim')]);

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
					self.colormap_matrices{j}{1} = feval(self.colormaps{j}{1},self.colormap_resolution);
					self.colormap_matrices{j}{2} = feval(self.colormaps{j}{2},self.colormap_resolution);
				else
					self.colormap_matrices{j} = feval(self.colormaps{j},self.colormap_resolution);
				end
			end
			self.refresh_slices;
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
		neg_idx = min(size(cmap{2},1),round((size(cmap{2},1)-1)*(-x-clim(1))/(clim(2)-clim(1)))+1);
		r = zeros(size(x));
		g = zeros(size(x));
		b = zeros(size(x));
		r(x>=0) = cmap{1}(idx(x>0),1);
		g(x>=0) = cmap{1}(idx(x>0),2);
		b(x>=0) = cmap{1}(idx(x>0),3);
		r(x<0) = cmap{2}(idx(x<0),1);
		g(x<0) = cmap{2}(idx(x<0),2);
		b(x<0) = cmap{2}(idx(x<0),3);
	else
		idx = min(size(cmap,1),round((size(cmap,1)-1)*(x-clim(1))/(clim(2)-clim(1)))+1);
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
