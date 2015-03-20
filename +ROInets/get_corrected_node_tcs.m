function [nodeData, voxelWeightings] = get_corrected_node_tcs(voxelData,    ...
                                                              spatialBasis, ...
                                                              protocol,     ...
                                                              timeCourseGenMethod)
%GET_CORRECTED_NODE_TCS correct ROI time-courses for source leakage
%
% methods for extracting node time courses using differing
% orthogonalisation protocols
%
% NODEDATA = GET_CORRECTED_NODE_TCS(VOXELDATA, PARCELFLAG, PROTOCOL, METHOD)
%   produces orthogonalised node time-courses NODEDATA from 
%   (nVoxels x nSamples) matrix of beamformed data VOXELDATA. VOXELDATA may
%   be passed in as an array, or as the filename of a saved .MAT file. 
%
%   PARCELFLAG is a logical array (nVoxels x nParcels) which identifies
%   the voxels in each parcel. Each column identifies the voxels 
%   making up each parcel. True entries indicate membership. 
%
%   PROTOCOL is a string to switch between various all-to-all 
%   orthogonalisation methods for source-spread correction. It can be:
%     'none'          - No correction applied. 
%     'symmetric'     - Apply orthogonalisation on the parcel time-courses.
%                       This produces orthonormal parcel time-courses
%                       which are as close as possible to the original
%                       time-courses.
%     'closest'       - Apply orthogonalisation on the parcel time-courses.
%                       Start as for the symmetric method, then converge to
%                       a (not orthonormal) orthogonal matrix which is as
%                       close as possible to the original time-courses. 
%     'householder'   - Orthogonalise using a more numerically stable
%                       alternative to the Gram-Schmidt process, dealing
%                       with ROI time-courses in a random order. 
%
%   METHOD chooses how the ROI time-course is created. It can be:
%     'PCA'   - take 1st PC of voxels 
%     'peakVoxel' - voxel with the maximum variance
%     note that the mean timecourse suffers issues relating to sign
%     ambiguities, and so is not currently an option.
%
% NODEDATA = GET_CORRECTED_NODE_TCS(VOXELDATA, SPATIALBASIS, PROTOCOL, METHOD)
%   it is also possible to use a spatial basis set (e.g. from group ICA) to
%   infer parcel time-courses. Each spatial map (held in columns) is a
%   whole-brain map - each map can be non-binary and maps may overlap.
%   SPATIALBASIS is therefore an (nVoxels x nSamples) matrix. 
%
%   Note that the spatial basis should be orthogonal, or nearly so. This is
%   because parcel data are computed one parcel at a time, without multiple
%   regression, which is valid for Data = SB * NodeData + e only if SB are
%   orthogonal. 
%
%   Options for METHOD are 
%      'spatialBasis' - The ROI time-course for each spatial map is the 1st 
%                       PC from all voxels, weighted by the spatial map. 
%
% [NODEDATA, VOXEL_WEIGHTINGS] = GET_CORRECTED_NODE_TCS(VOXELDATA, SPATIALBASIS, PROTOCOL, METHOD)
%    will return the relative weighting of voxels over the brain, used to
%    construct the time-course for each ROI. 


% Spatial basis algorithm:
% 1. Scale all maps so they have a positive peak of height 1
% 
% For each map in the spatial basis:
% 
%     2. Variance-normalise the clean voxel time-series. (This is now turned off.) 
%        Weight by the clean spatial map. 
% 
%     3. Perform a PCA. Extract the coefficients of the first PC to represent the node time-course
% 
%     4. Re-weight the node time-course to match the variance in the ROI. 
    

%	Copyright 2014 OHBA, FMRIB
%	This program is free software: you can redirstribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 375 $
%	$LastChangedDate: 2015-01-12 18:21:57 +0000 (Mon, 12 Jan 2015) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough and Stephen Smith.

%% find representative time courses for each parcel
% include quick hack to allow saved matrices to be passed in, to save on
% memory
if ischar(voxelData) && strcmp(voxelData(end-3:end), '.mat') && exist(voxelData, 'file'), 
    tmp = load(voxelData);
    ff = fieldnames(tmp);
    voxelData = tmp.(ff{1});
    clear tmp;
elseif ischar(voxelData),
    error([mfilename ':FileInputError'],               ...
          ['Did not recognise the filename input. \n', ...
           '  Ensure the file exists and has a .mat extension. ']);
else
    % carry on
end%if

nParcels = ROInets.cols(spatialBasis);

ft_progress('init', 'text', '');

switch lower(timeCourseGenMethod)
    case 'pca'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'It will be binarised. \n']);
        end%if
        spatialBasis = logical(spatialBasis);    
        % check that each voxel is only a member of one parcel
        assert(~any(ROInets.row_sum(spatialBasis) > 1), ...
               [mfilename ':MultipleParcelOccupancy'], ...
               'Each voxel can be a member of at most one parcel. \n');

        % demean and variance normalise each voxel - 20 May 2014 Don't
        % variance normalise as MEG data are very smooth, and power is a
        % good indication of true sources
        temporalSTD     = max(std(voxelData, [], 2), eps);
        voxelDataScaled = ROInets.demean(voxelData, 2);
        clear voxelData
        
        % pre-allocate PCA weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
        
        % perform PCA on each parcel and select 1st PC scores to represent
        % parcel
        for iParcel = nParcels:-1:1,
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding PCA time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
                
            thisMask = spatialBasis(:, iParcel);
            if any(thisMask), % non-zero
                parcelData = voxelDataScaled(thisMask, :);
                
                [U, S, V]  = ROInets.fast_svds(parcelData, 1);
                PCAscores  = S * V';
                
                % restore sign and scaling of parcel time-series
                % U indicates the weight with which each voxel in the
                % parcel contributes to the 1st PC
                TSsign          = sign(mean(U));
                relVoxelWeights = abs(U) ./ sum(abs(U)); % normalise the linear combination
                % weight the temporal STDs from the ROI by the proportion used in 1st PC
                TSscale         = dot(relVoxelWeights, temporalSTD(thisMask)); 
                nodeTS          = TSsign .*                               ...
                                  (TSscale / max(std(PCAscores), eps)) .* ... 
                                  PCAscores;
                      
                % return the linear operator which is applied to the data
                % to retrieve the nodeTS
                voxelWeightings(thisMask, iParcel) = TSsign .* ...
                                                     (TSscale / max(std(PCAscores), eps)) ...
                                                     .* U';
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(voxelDataScaled));
            end%if
            
            nodeDataOrig(iParcel,:) = nodeTS;
        end%for
        
        clear parcelData voxelDataScaled

    case 'peakvoxel'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'It will be binarised. \n']);
        end%if
        spatialBasis = logical(spatialBasis);
        
        % find rms power in each voxel
        voxelPower = sqrt(ROInets.row_sum(voxelData.^2) ./ ...
                          ROInets.cols(voxelData));
                      
        % pre-allocate weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
                      
        % take peak voxel in each parcel
        for iParcel = nParcels:-1:1,
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding peak voxel time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            thisMask = spatialBasis(:, iParcel);
            
            if any(thisMask), % non-zero
                % find index of voxel with max power
                thisParcPower            = voxelPower;
                thisParcPower(~thisMask) = 0;
                [~, maxPowerInd]         = max(thisParcPower);
                
                % select voxel timecourse
                nodeDataOrig(iParcel,:) = voxelData(maxPowerInd,:);
                
                % save which voxel was used
                voxelWeightings(maxPowerInd, iParcel) = 1;
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeDataOrig(iParcel,:) = zeros(1, ROInets.cols(voxelData));
            end%if
        end%loop over parcels
        
        clear voxelData parcelData
        
    case 'spatialbasis'
        % scale group maps so all have a positive peak of height 1
        % in case there is a very noisy outlier, choose the sign from the
        % top 5% of magnitudes
        top5pcInd = abs(spatialBasis) >=                        ...
                         repmat(prctile(abs(spatialBasis), 95), ...
                                [ROInets.rows(spatialBasis), 1]);
        for iParcel = nParcels:-1:1,
            mapSign(iParcel) = sign(mean(...
                              spatialBasis(top5pcInd(:,iParcel), iParcel)));
        end%for
        scaledSpatialMaps = ROInets.scale_cols(spatialBasis, ...
                                   mapSign ./                ...
                                   max(max(abs(spatialBasis), [], 1), eps));
        
        % estimate temporal-STD for normalisation
        temporalSTD = max(std(voxelData, [], 2), eps);
                               
        % find time-course for each spatial basis map
        for iParcel = nParcels:-1:1, % allocate memory on the fly
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding spatial basis time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            % extract the spatial map of interest
            thisMap     = scaledSpatialMaps(:, iParcel);
            parcelMask  = logical(thisMap);
                        
            % variance-normalise all voxels to remove influence of
            % outliers. - remove this step 20 May 2014 for MEG as data are
            % smooth and little risk of high-power outliers. Also, power is
            % a good indicator of sensible signal. 
            % Weight all voxels by the spatial map in question
            weightedTS  = ROInets.scale_rows(voxelData, thisMap);
            
            % perform svd and take scores of 1st PC as the node time-series
            % U is nVoxels by nComponents - the basis transformation
            % S*V holds nComponents by time sets of PCA scores - the 
            % timeseries data in the new basis
            [U, S, V]   = ROInets.fast_svds(weightedTS(parcelMask,:), 1);
            clear weightedTS
            
            PCAscores   = S * V';
            maskThresh  = 0.5; % 0.5 is a decent arbitrary threshold chosen by Steve Smith and MJ after playing with various maps.
            thisMask    = thisMap(parcelMask) > maskThresh;   
            
            if any(thisMask), % the mask is non-zero
                % U is the basis by which voxels in the mask are weighted
                % to form the scores of the 1st PC
                relativeWeighting = abs(U(thisMask)) ./ ...
                                    sum(abs(U(thisMask)));
                
                TSsign  = sign(mean(U(thisMask)));
                TSscale = dot(relativeWeighting, temporalSTD(thisMask));       
                nodeTS  = TSsign .*                               ...
                          (TSscale / max(std(PCAscores), eps)) .* ...      
                          PCAscores;
                      
                % for Mark: this is the linear operator which is applied to
                % the voxel data to get nodeTS.
                voxelWeightings(parcelMask,iParcel) = TSsign .* ...
                                             (TSscale / max(std(PCAscores), eps)) ...
                                             .* (U' .* thisMap(parcelMask)');
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(weightedTS));
                voxelWeightings(~thisMask, iParcel) = zeros(length(thisMask), 1);
            end%if
            
            nodeDataOrig(iParcel, :) = nodeTS;
            
        end%loop over parcels
        
        clear voxelData 
        
        
    otherwise
        error([mfilename ':UnrecognisedTimeCourseMethod'],            ...
              ['Unrecognised method for finding ROI time-course. \n', ...
               'Expected ''PCA'', ''spatialBasis'', or ''mean''. \n']);
end%switch

ft_progress('close');

%% Do orthogonalisation
nodeData = ROInets.remove_source_leakage(nodeDataOrig, protocol);

end%get_orthogonalised_node_tcs
% [EOF]