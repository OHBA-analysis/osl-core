classdef parcellation < handle
	% Class-based representation of parcellations
	%
	% Designed to work with parcellations on OSL standard masks only
	% e.g. currently doesn't work with HCP/fieldtrip coordinates
	% 
	% KEY REPRESENTATIONS
	%  - p.weight_mask - XYZ x Parcels
	%
	% KEY MANIPULATIONS
	% - p.to_vol - Convert a matrix to a 4D volume
	% - p.to_matrix - Convert a 4D volume to a matrix
	%
	% Try p.plot or p.fslview to display regions

	properties
	    weight_mask % XYZ x Parcels matrix of parcel weights
	    template_mask % Values for each voxel in the template
	    template_coordinates % MNI coordinates for each point in the mask
	    template_fname % Filename of standard template, useful for copying qform/sform matrix 
	    labels % Names of each ROI in the parcellation
	end

	% These properties are dependent on the weight mask, so are only ever set indirectly
	% They are protected rather than dependent so they are only computed once
	properties(SetAccess = protected)
		is_weighted % 'true' if this is a weighted parcellation
		is_overlapping % 'true' if any voxels belong to more than one parcel
		resolution % Mask resolution (mm)
	end

	properties(Dependent)
		n_parcels % Number of parcels
		n_voxels % Number of voxels
	end

	methods

		function self = parcellation(input_mask,labels,template)
			% Take in a mask or a file storing a mask, convert it to a Parcellation
			%
			% INPUTS
			%	- input_mask specifying the parcels themselves
			%	- labels specifying the ROI names (one for each parcel)
			% 	- template specifying the voxel coordinates (currently custom templates not supported, only parcellations based on OSL standard masks)
			%
			% INPUT_MASK formats that can be loaded
			% - XYZ x Parcels 
			% - Vox x Parcels 
			% - Vox x 1 
			% - XYZ x 1
			% - nii file name (nifti file with one of the above formats)
			% - mat file name (containing variables 'mask' for input mask and 'long names' for labels)
			% See set.weight_mask for how input matrix is reshaped and analyzed
			%
			% LABELS formats that can be loaded
			% - Cell array of strings
			% - mat file name containing variable 'labels'
			% - txt file with each parcel name on a single row
			%
			% If no labels are specified, and the input file is not a mat file, default labels will be generated either by
			% - If a .nii file was loaded, looking for a .txt file with the same file name
			% - Otherwise, placeholder labels 'ROI 1', 'ROI 2' etc. will be generated
					
			% If file is not present, try looking in the OSL parcellations folder
			if ischar(input_mask) && ~exist(input_mask)
				input_mask = fullfile(getenv('OSLDIR'),'parcellations',input_mask);
			end

			if nargin < 2 || isempty(labels) 
				labels = [];
			else
				if ~exist(labels)
					labels = fullfile(getenv('OSLDIR'),'parcellations',labels);
				end

				[~,~,ext] = fileparts(labels);

				switch ext
					case '.txt'
						labels = importdata(labels);
					case '.mat'
						d = load(labels);
						labels = d.labels;
					otherwise
						error('OSL:UnrecognizedExtension',sprintf('Unrecognized label extension ''%s'', must end with .txt or .mat',ext));
				end
			end

			% Load the input file
			if ischar(input_mask) && ~isempty(strfind(input_mask,'nii'))  % nifti input
				
				% See if there is a corresponding label file in the same directory
				% Do this first because input_mask gets converted from filename to mask data
				if isempty(labels)
					label_fname = strrep(input_mask,'nii','txt');
					if exist(label_fname)
						labels = importdata(label_fname);
					end
				end

				input_mask = readnii(input_mask);

			elseif ischar(input_mask) % matlab input
				d = load(input_mask);
				input_mask = d.mask;

				% Support legacy files that store parcellations as cell array of ROI coordinates
				if isfield(d,'roi_coords') 
					input_mask = zeros(size(d.mask,1),1);
					for j = 1:length(d.roi_coords)
						input_mask(ismember(d.mask,d.roi_coords{j},'rows')) = j;
					end
				end

				% If no label file is provided, try and get the labels from the .mat file
				if isempty(labels)
					if isfield(d,'labels')
						labels = d.labels;
					elseif isfield(d,'long_names')
						labels = d.long_names;
					end
				end
			elseif length(input_mask) == 1 % Can enter a spatial resolution to retrieve the whole brain mask i.e. 1 parcel
				[~,self.template_fname,self.template_mask] = self.guess_template(input_mask);
				input_mask = ones(sum(self.template_mask(:)>0),1);
			end

			% Guess and load template if required
			if nargin < 3 || isempty(template) 
				[self.resolution,self.template_fname,self.template_mask] = self.guess_template(input_mask);
				self.template_coordinates = osl_mnimask2mnicoords(self.template_fname);
			else
				error('Specifying a nonstandard template not yet supported')
			end
			
			if ndims(self.template_mask) < 3
				[~,~,mask_1] = self.guess_template(input_mask); % If user gave their own template on a standard grid size
				self.template_mask = matrix2vols(self.template_mask,mask_1);
			end

			self.weight_mask = input_mask;

			% Assign default labels if required
			if isempty(labels)
				labels = cell(self.n_parcels,1);
				for j = 1:self.n_parcels
					labels{j} = sprintf('ROI %d',j);
				end
			else
				assert(iscell(labels),'Manually specified labels must be a cell array of ROI names or a file name');
				assert(isvector(labels),'Labels must be a cell vector of strings, not a matrix');
				assert(length(labels)==self.n_parcels,sprintf('Must have one manually specified label for each parcel (%d provided, %d required)',length(labels),self.n_parcels));
			end

			self.labels = labels(:);
		end

		function set.weight_mask(self,mask)
			% Possible mask sizes
			% - Vox x 1 
			% - Vox x Parcels 
			% - XYZ x 1 (3D matrix, unweighted non-overlapping parcel if integers OR single weighted non-overlapping parcel if values are not integers)
			% - XYZ x Parcels (4D matrix, unweighted if all values 0 or 1, overlapping if multiple assingmnet)

			% First, convert Vox x 1 or Vox x Parcels to XYZ x 1 or XYZ x Parcels

			if ndims(mask) < 3
				mask = matrix2vols(mask,self.template_mask);
			end
			
			if ndims(mask) == 3
				all_integers = all(mod(self.weight_mask(:),1)==0); % Do we only have integers?
				vals = unique(self.weight_mask(:));
				self.is_weighted = ~(all_integers && all(diff(vals) == 1)); % If we read in 3 dimensions, continuous integers mean interpreted parcel index rather than weight
				
				if ~self.is_weighted
					mask = self.integers_to_masks(mask);
				end
			elseif ndims(mask) == 4
				self.is_weighted = ~all(ismember(mask(:),[0 1])); % If we read in 4 dimensions, then any non-binary value means its weighted
			end

			% Finally, check if overlapping
			parcels_per_voxel = sum(mask~=0,4);
			if max(parcels_per_voxel(:)) > 1
				self.is_overlapping = true;
			else
				self.is_overlapping = false;
			end

			self.weight_mask = mask;
		end

	    function parcelflag = parcelflag(self,binarize)
			% Return the parcelflags for ROInets.get_node_tcs
			if nargin < 2 || isempty(binarize) 
				binarize = false;
			end
			if binarize
				p = self.binarize;
			else
				p = self.weight_mask;
			end

			parcelflag = self.to_matrix(p);
		end

		function dat4 = to_vol(self,dat2);
			% Convert Vox x Frames or Parcels x Frames, to XYZ x Frames
			%
			% Commonly used to produce a matrix suitable for use with fslview

			% Ensure matrix is in correct orientation
			% If matrix is 
			if size(dat2,1) ~= self.n_voxels && size(dat2,1) ~= self.n_parcels % If the first dimension is neither n_voxels nor n_parcels
				dat2 = dat2.';
			end

			if size(dat2,1) == self.n_parcels
				% This means we need to expand to voxels first
				d2 = nan(self.n_voxels,size(dat2,2));
				m = self.value_vector; % Map voxels to parcels
				for k = 1:size(dat2,2)
					for j = 1:self.n_parcels
						d2(m==j,k) = dat2(j,k);
					end
				end
				dat2 = d2;
			elseif ~size(dat2,1) == self.n_voxels
				error('Unsupported dimension')
			end

			dat4 = matrix2vols(dat2,self.template_mask);
		end

		function dat2 = to_matrix(self,dat4)
			% Convert XYZ x Frames to Vox x Frames
			% If no matrix is provided, it will return the matrix
			% representation of the parcellation
			if nargin < 2 || isempty(dat4) 
				dat4 = self.weight_matrix;
			end
			
			dat2 = vols2matrix(dat4,self.template_mask);
		end

		function p = value_vector(self)
			% Return a Vox x 1 vector where value is parcel index (binarizes if necessary)
			if self.is_weighted || self.is_overlapping
				fprintf(2,'Warning - parcellation is being binarized\n')
			end
			p = self.to_matrix(self.binarize);
            p = bsxfun(@times,p,1:self.n_parcels);
			p = sum(p,2);
		end

		function weight_mask = remove_weights(self)
			% Return a weight mask with weightings removed, but overlap still permitted
			weight_mask = +(self.weight_mask ~= 0);
		end
		
		function weight_mask = remove_overlap(self)
			% Return a weight mask with no overlap, but original weights preserved
			% e.g. If voxel 1 belongs to parcel 1 (0.5) and parcel 2 (0.25) and it's assigned
			% to parcel 1, then the new weights will be parcel 1 (0.5) and parcel 2 (0)
			weight_mask = self.weight_mask.*self.binarize();
		end

		function weight_mask = binarize(self)
			% Return a weight mask with no weights and no overlap
			% binarize() is equivalent to remove_overlap() followed by remove_weights()
			% WARNING - it is NOT equivalent to remove_weights() followed by remove_overlap()
			p = self.to_matrix(self.weight_mask);
			[~,assignments] = max(p,[],2);
			assignments(all(p==0,2)) = 0; % Voxels x 1 with value indicating parcel
			assignment_vol = self.to_vol(assignments); % XYZ x 1 with values indicating parcel
			weight_mask = self.integers_to_masks(assignment_vol); % XYZ x n_parcels
		end

		% PROPERTIES
		function n = get.n_parcels(self)
			n = size(self.weight_mask,4);
		end

		function n = get.n_voxels(self)
			n = size(self.template_coordinates,1);
		end

		function plot(self)
			self.show_parcellation();
		end

		function fslview(self,single_volume)
			% Display this parcellation using fslview
			% If single_volume == true, then there will be one binarized
			% volume where the value indicates the parcel
			if nargin < 2 || isempty(single_volume) 
				single_volume = false;
			end
			
			% Display the parcellation using fslview
			if single_volume
				p = self.to_vol(self.value_vector);
				clim = [0 self.n_parcels];
			else
				p = self.weight_mask;
				clim = [0 max(self.weight_mask(:))];
			end

			fname = self.savenii(p);
			fslview(fname,clim);
			delete(fname)
		end

		function fname = savenii(self,data,fname)
			% Save a nii file, with qform/xform copied from the original mask file
			if nargin < 3 || isempty(fname) 
				fname = tempname('.');
			end
					
		    [~,~,scales] = read_avw(self.template_fname);
		    save_avw(data,fname,'f',scales);
		    system(['fslcpgeom ' self.template_fname ' ' fname ' -d']);
		end

		function [coords,weights] = roi_coords(self)
			% Retrive the voxel coordinates associated with each ROI
			%
			% Returns a cell array where the cell index is the ROI
			% and the contents is a list of voxel coordinates associated
			% with that ROI

			w = self.to_matrix(self.weight_mask);
			for j = 1:self.n_parcels
				idx = w(:,j) ~= 0;
				coords{j} = self.template_coordinates(idx,:);
				weights{j} = w(idx,j);
			end

		end

	end

	methods (Static)
		function d4 = integers_to_masks(d3)
			% Input - d3 = XYZ x 1 where all values are integers
			% Output - d4 = XYZ x n_parcels (binary parcellation)
			assert(ndims(d3) == 3,'Input must be XYZ x 1')
			assert(all(mod(d3(:),1)==0),'Input must only contain integers'); 
			n_parcels = max(d3(:));
			d4 = zeros([size(d3) n_parcels]);
			for j = 1:n_parcels
				d4(:,:,:,j) = d3 == j;
			end
		end

		function [spatial_res,mask_fname,mask] = guess_template(m)
			% Given a matrix of parcellations or similar, guess and return the matching mask
			% m could be
			% - The spatial resolution desired (e.g. m=8)
			% - A column vector of parcels (voxels x 1)
			% - A matrix of parcels (voxels x parcels)
			% - A volume (XYZ)
			% - A volume-based parcellation (XYZ x parcels)
			% 
			% Return values
			% mask - The result of reading the mask file using readnii()
			% mask_fname - The filename for the standard mask
			% spatial_res - The guessed spatial resolution

	    
			if isa(m,'meeg')
				m = m(:,1,1);
			end 
			% Given a 4D matrix, it's probably XYZ x parcels, so just keep the first volume
			if ndims(m) == 4
				m = squeeze(m(:,:,:,1));
			end

			% If its a 2D matrix, its probably voxels x parcels, so just keep the first column
			if ndims(m) == 2
				m = m(:,1);
			end

			mask_res = [1:15];

			% Commands below generate the sizes listed in this file
			% OSLDIR = getenv('OSLDIR');
			% mask_dim = [];
			% mask_vox = [];
			% for j = 1:length(mask_res)
			% 	a=readnii(sprintf([OSLDIR '/std_masks/MNI152_T1_%dmm_brain.nii.gz'],mask_res(j)));
			% 	mask_dim(j,:) = size(a);
			% 	mask_vox(j) = sum(a(:)~=0);
			% end

			mask_dim = [...
			   182   218   182;...
			    91   109    91;...
			    61    73    61;...
			    46    55    46;...
			    36    44    36;...
			    30    36    30;...
			    26    31    26;...
			    23    27    23;...
			    20    24    20;...
			    18    22    18;...
			    17    20    17;...
			    15    18    15;...
			    14    17    14;...
			    13    16    13;...
			    12    15    12];

			mask_vox = [1827095,228453,67693,28549,14641,8471,5339,3559,2518,1821,1379,1065,833,676,544];

			if numel(m) == 1
				idx = find(mask_res == m);
			elseif size(m,2) == 1 % If column vector given
				idx = find(mask_vox == size(m,1));
				if isempty(idx)
					idx = find(prod(mask_dim,2) == size(m,1));
					if isempty(idx)
						error('No masks have %d voxels',size(m,1));
					end
				end
			else
				idx = find(all(bsxfun(@eq,mask_dim,size(m)),2));
				if isempty(idx)
					size(m)			
					error('No masks have the right dimension')
				end
			end

			spatial_res = mask_res(idx);

			if nargout > 1
				OSLDIR = getenv('OSLDIR');
				mask_fname = [OSLDIR '/std_masks/MNI152_T1_' num2str(spatial_res) 'mm_brain.nii.gz'];
			end

			if nargout > 2
				mask = readnii(mask_fname);
			end
		end

	end

end
