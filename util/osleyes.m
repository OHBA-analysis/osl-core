function osleyes(niifiles)
	% To behave properly - we have a single MNI coordinate system
	% which is where the crosshairs are pointed

	% For consistency - *always* flip left-right

	niifiles = {'std_masks/MNI152_T1_8mm_brain.nii.gz'}

	% Set up three axes
	h.fig = figure;
	h.ax(1) = axes;
	h.ax(2) = axes;
	h.ax(3) = axes;
	set([h.fig h.ax],'Units','normalized');

	for j = 1:length(niifiles)
		imgfile = load_untouch_nii(niifiles{j}); % Include applying xform/qform
		h.imgdata{j} = imgfile.img(:,:,:); % Flip left to right
		h.coord{j} = get_coords(imgfile.hdr);
		h.img{j} = render_image(h.imgdata{j},h.coord{j},h.ax);
	end

	%hold(ax,'on');
	% h_img(1) = imagesc(permute(img{1}(10,:,:),[2 3 1]),'Parent',ax(1),'HitTest','off')
	% h_img(2) = imagesc(permute(img{1}(:,10,:),[1 3 2]),'Parent',ax(2),'HitTest','off')
	% h_img(3) = imagesc(permute(img{1}(:,:,10),[1 2 3]),'Parent',ax(3),'HitTest','off')

	set(h.ax(1),'Position',[0 0.2 0.3 0.8],'View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal');
	set(h.ax(2),'Position',[0.35 0.2 0.3 0.8],'View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
	set(h.ax(3),'Position',[0.7 0.2 0.3 0.8],'View', [0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');

	hold(h.ax(1),'on');
	hold(h.ax(2),'on');
	hold(h.ax(3),'on');

	h.h_crosshair(1) = plot(h.ax(1),NaN,NaN,'g','HitTest','off');
	h.h_crosshair(2) = plot(h.ax(2),NaN,NaN,'g','HitTest','off');
	h.h_crosshair(3) = plot(h.ax(3),NaN,NaN,'g','HitTest','off');
	set(h.ax,'ButtonDownFcn',@axis_click)
	colormap bone

	h.current_point = [10 10 10];

	guidata(h.fig,h)

	update_slices(h.fig)

end

function h_img = render_image(img,coord,ax)
	h_img(1) = imagesc(coord.y,coord.z,permute(img(1,:,:),[2 3 1])','Parent',ax(1),'HitTest','off');
	h_img(2) = imagesc(coord.x,coord.z,permute(img(:,1,:),[1 3 2])','Parent',ax(2),'HitTest','off');
	h_img(3) = imagesc(coord.x,coord.y,permute(img(:,:,1),[1 2 3])','Parent',ax(3),'HitTest','off');
end

function axis_click(a,b)
	h = guidata(a);
	p = get(a,'CurrentPoint');

	if a == h.ax(1)
		h.current_point(2:3) = p(1,1:2);
	elseif a == h.ax(2)
		h.current_point([1,3]) = p(1,1:2);
	else
		h.current_point(1:2) = p(1,1:2);
	end
	h.current_point
	guidata(h.fig,h);
	update_slices(h.fig)
end

function update_slices(fig)
	% Given crosshair position, update the slices
	h = guidata(fig);
	p = h.current_point; % Current point in 3D

	set(h.h_crosshair(1),'XData',[get(h.ax(1),'XLim') NaN p(2) p(2)],'YData',[p(3) p(3) NaN get(h.ax(1),'YLim')]);
	set(h.h_crosshair(2),'XData',[get(h.ax(2),'XLim') NaN p(1) p(1)],'YData',[p(3) p(3) NaN get(h.ax(2),'YLim')]);
	set(h.h_crosshair(3),'XData',[get(h.ax(3),'XLim') NaN p(1) p(1)],'YData',[p(2) p(2) NaN get(h.ax(3),'YLim')]);

	% Now update each slice

	for j = 1:length(h.img)
		[~,idx(1)] = min(abs(h.coord{j}.x-p(1)));
		[~,idx(2)] = min(abs(h.coord{j}.y-p(2)));
		[~,idx(3)] = min(abs(h.coord{j}.z-p(3)));
		set(h.img{j}(1),'CData',permute(h.imgdata{j}(idx(1),:,:),[2 3 1])');
		set(h.img{j}(2),'CData',permute(h.imgdata{j}(:,idx(2),:),[1 3 2])');
		set(h.img{j}(3),'CData',permute(h.imgdata{j}(:,:,idx(3)),[1 2 3])');
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

