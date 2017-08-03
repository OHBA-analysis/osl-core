classdef osleyes < handle

	properties
		clims = {}; % Values outside range are transparent!
		colormaps= {};
		current_point = [1 1 1]
	end

	properties(SetAccess=protected)
		fig
		ax
		nii
		coord
		h_img = {} % Cell array of img handles associated with each nii file (there are 3)
		h_crosshair
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
				self.coord{j} = get_coords(self.nii{j}.hdr);
				self.h_img{j}(1) = imagesc(self.coord{j}.y,self.coord{j}.z,permute(self.nii{j}.img(1,:,:),[2 3 1])','Parent',self.ax(1),'HitTest','off');
				self.h_img{j}(2) = imagesc(self.coord{j}.x,self.coord{j}.z,permute(self.nii{j}.img(:,1,:),[1 3 2])','Parent',self.ax(2),'HitTest','off');
				self.h_img{j}(3) = imagesc(self.coord{j}.x,self.coord{j}.y,permute(self.nii{j}.img(:,:,1),[1 2 3])','Parent',self.ax(3),'HitTest','off');
				self.clims{j} = [min(self.nii{j}.img(:)),max(self.nii{j}.img(:))];
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

			colormap(self.ax(1),'bone')
			colormap(self.ax(2),'bone')
			colormap(self.ax(3),'bone')

			self.refresh_slices()

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

				d1 = permute(self.nii{j}.img(idx(1),:,:),[2 3 1])';
				d2 = permute(self.nii{j}.img(:,idx(2),:),[1 3 2])';
				d3 = permute(self.nii{j}.img(:,:,idx(3)),[1 2 3])';
				set(self.h_img{j}(1),'CData',d1,'AlphaData',+(d1>self.clims{j}(1)));
				set(self.h_img{j}(2),'CData',d2,'AlphaData',+(d2>self.clims{j}(1)));
				set(self.h_img{j}(3),'CData',d3,'AlphaData',+(d3>self.clims{j}(1)));
			end
		end

		function refresh_colors(self)

		end

	end

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

