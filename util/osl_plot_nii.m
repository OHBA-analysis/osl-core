function h = osl_plot_nii(fnames,thresholds,colormaps,anatomical)
	% Display a NIFTI file on top of an (optionally) automatically selected 
	% anatomical image
	%
	% This file automatically selects an anatomical image, and displays 
	% the input nii's as overlays on this image. An osleyes() object will be 
	% instantiated and returned by this function.
	%
	% INPUTS
	% - fnames : Name of an image, or a cell array with m images
	% - Thresholds : Optional m x 2 matrix or cell array of length 2 vectors with thresholds
	% - colormaps : Optional m x 1 cell array of colormaps to pass to osleyes
	% - anatomical : File name of a base image to display as the bottom later
	%
	% Note that the main intended usage of osl_plot_nii as opposed to osleyes() is
	% to automatically select and include the anatomical image. It is generally 
	% recommended that properties like colormaps and color ranges are set 
	% by interacting with the osleyes object directly

    if ~iscell(fnames)
		fnames = {fnames};
    end
    
	if nargin < 4 || isempty(anatomical) 
		vol = osl_load_nii(fnames{1});
		[~,anatomical,~] = parcellation.guess_template(vol);
	end

	if nargin < 3 || isempty(colormaps) 
		colormaps = [];
	end
	
	if nargin < 2 || isempty(thresholds) 
		thresholds = [];
	else
		if isnumeric(thresholds)
			tmp = thresholds;
			thresholds = {};
			for j = 1:size(tmp)
				thresholds{j} = tmp(j,:);
			end
		end
    end

	fnames = [{anatomical},fnames];

	h = osleyes(fnames);

	% Now set the thresholds and colormaps

	if ~isempty(thresholds)
		h.clims(2:end) = thresholds;
	end

	if ~isempty(colormaps)
		h.colormaps(2:end) = colormaps;
	end
