function show_parcellation(p)
	% Interactive tool to view ROIs given a p and corresponding mask
	
	plot_parcel_volume = false;

	mask = p.template_mask;
	rois.masks = p.weight_mask;
	rois.long_names = p.labels;
	mni_coords = p.template_coordinates;
 %    mask       = load_untouch_nii(mask);
 %    mni_coords = nifti_mask2coord(mask);
	% rois       = load(parcel);


	f=figure;
	shp = alphaShape(mni_coords(:,1),mni_coords(:,2),mni_coords(:,3));
	h_brain = shp.plot;
	h_brain.FaceAlpha = 0.1;
	h_brain.FaceColor = [0.5 0.5 1];
	h_brain.EdgeColor = 'w';
	h_brain.EdgeAlpha = 0.2;
	ax = gca;

	hold on

	nstr = {};
	for j = 1:length(rois.long_names)
		nstr{j} = sprintf('%d - %s',j,rois.long_names{j});
	end

	%a = uicontrol(f,'Style','popupmenu','String',arrayfun(@(x) sprintf('%d',x),1:size(rois.masks,4),'UniformOutput',false),'Value',1);
	a = uicontrol(f,'Style','popupmenu','String',nstr,'Value',1);

	a.Position(3) = 100;

	h_roi = scatter3(NaN,NaN,NaN,30,'r');
	h_roi_vol = [];

	a.Callback = @(a,b,c) draw_roi(rois,a.Value,h_roi);
	a.Callback(a); % Show the first ROI immediately


	function draw_roi(rois,idx,h_roi)
	    roi = rois.masks(:,:,:,idx); % 0 and 1 for grid values
	    roi = logical(vols2matrix(roi,p.template_mask));
	    % roi = roi(logical(mask.img(:)));
	    h_roi.XData = mni_coords(roi,1);
	    h_roi.YData = mni_coords(roi,2);
	    h_roi.ZData = mni_coords(roi,3);

	    if plot_parcel_volume
		    delete(h_roi_vol)
		    shp = alphaShape(mni_coords(roi,1),mni_coords(roi,2),mni_coords(roi,3));
		    h_roi_vol = shp.plot;
		    h_roi_vol.FaceAlpha = 0.3;
		    h_roi_vol.FaceColor = [1 0 0];
		    h_roi_vol.LineStyle = 'none';
		end
		    
	end
end
