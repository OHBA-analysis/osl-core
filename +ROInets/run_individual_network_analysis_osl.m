function CorrMats = run_individual_network_analysis_osl(oat,               ...
                                                        ReconResultsfName, ...
                                                        spatialBasis,      ...
                                                        Settings,          ...
                                                        iSession,          ...
                                                        resultsSaveName)
%RUN_INDIVIDUAL_CORRELATION_ANALYSIS_OSL calleld by osl_network_analysis
%
% CORRELATION_MATS = run_individual_correlation_analysis_osl(OAT,
%  RECON_RESULTS_FILE_NAME, SPATIAL_BASIS, SETTINGS, SESSION_NAME, RESULTS_SAVE_NAME) performs
%   a network correlation analysis on the session with source recon results
%   file RECON_RESULTS_FILE_NAME, from a source-recon run on OAT. ROIs are
%   defined by a spatial basis set or set of binary ROI masks
%   SPATIAL_BASIS, with settings (oil.ROInetworks) held in SETTINGS. (Type
%   `help osl_network_analysis' for more info.) SESSION_NAME defines the
%   name of the saved node data files (the corected node time-courses).
%   RESULTS_SAVE_NAME sets the location where the individual session
%   correlation matrices will be saved to disc. 


%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
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
%	$Revision: 263 $
%	$LastChangedDate: 2014-10-23 11:30:39 +0100 (Thu, 23 Oct 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 23:39:39

%% Unpack settings
sessionName     = Settings.sessionNames{iSession};
doRegularize    = Settings.Regularize.do;
nH0Iter         = Settings.nEmpiricalSamples;
ARorder         = Settings.ARmodelOrder;
tcsMethod       = Settings.timecourseCreationMethod;

switch lower(Settings.leakageCorrectionMethod),
    case 'pairwise',
        % this is a bit of a fudge at the moment. We perform the
        % orthogonalisation and correlation steps together
        protocol = 'none';
        
    otherwise
        protocol = Settings.leakageCorrectionMethod;
end%switch

%% Extracting timecourses

[~, ReconResultsLoadfName]        = fileparts(ReconResultsfName);
[voxelDataBroadBand, allTime, Fs] = MEGsim.get_source_data(oat, ReconResultsLoadfName);

% select time range from user-defined input
[time, timeInd, timeRange] = time_range(allTime, Settings.timeRange, iSession);
% apply to source data, too.
voxelDataBroadBand(:, ~timeInd) = [];
nTimeSamples                    = ROInets.cols(voxelDataBroadBand);
nNodes                          = ROInets.cols(spatialBasis);

assert(sum(timeInd) == nTimeSamples, ...
       [mfilename ':InconsistentTimeIndices'], ...
       'There''s been an internal error in selecting the time range. \n');
   
%% Apply ROIs to save RAM
assert(isequal(ROInets.rows(spatialBasis), ROInets.rows(voxelDataBroadBand)), ...
       [mfilename ':ROIDataMisMatch'],                                        ...
       'Number of voxels in data and ROI basis set do not match. \n');
   
% check we have a binary spatial basis
isBinaryMask = all(0 == spatialBasis(:) | 1 == spatialBasis(:));
if isBinaryMask,
    fprintf('Applying ROIs. \n');
    allROImask = logical(ROInets.row_sum(spatialBasis));
    % remove irrelevant voxels
    voxelDataBroadBand(~allROImask, :) = [];
    spatialBasis(~allROImask, :)       = [];
else
    allROImask = true(ROInets.rows(spatialBasis), 1);
end%if

%% Fit an AR model for the construction of empirical H0 data
if nH0Iter,
    [ARmodel.coeffs,           ...
     ARmodel.varianceEstimate, ...
     ARmodel.partial_coeffs  ] = ROInets.estimate_AR_coeffs(voxelDataBroadBand, ...
                                                            ARorder);
else
    ARmodel = [];
end%if

%% Perform separate analyses in each frequency band
for iFreq = Settings.nFreqBands:-1:1,
    fprintf(' Correlation analysis for frequency band %d of %d. \n', ...
            iFreq, Settings.nFreqBands);
        
    %% bandpass filter the data
    voxelData = zeros(size(voxelDataBroadBand)); % pre-allocate memory for nested functions
    [Filter, bandName] = filter_voxel_data(Settings, Fs, iFreq);
    
    %% find node time-courses, 
    % applying symmetric orthogonalisation if necessary
    fprintf([' Finding time courses with method %s;\n', ...
             ' Removing source leakage with orthogonalisation method %s.\n'], ...
            tcsMethod, Settings.leakageCorrectionMethod);
    

    [nodeData, voxelWeightings] = get_corrected_node_tcs(spatialBasis, ...
                                                         protocol,     ...
                                                         tcsMethod);
    
    save_corrected_timecourse_results(nodeData, allROImask, ...
                                      voxelWeightings, Settings,          ...
                                      sessionName, protocol, bandName);
      
    %% correlation analysis
    if strcmpi(Settings.leakageCorrectionMethod, 'pairwise'),
        CorrMats{iFreq} = ROInets.do_pairwise_calculation(nodeData, ...
                                                          time,     ...
                                                          Settings.EnvelopeParams);
        % free up some memory
        clear nodeData
        
    else
        fprintf('  Enveloping \n');
        
        % take power envelopes
        nodeEnv = ROInets.envelope_data(nodeData,        ...
                                        time,            ...
                                        Settings.EnvelopeParams);
        
        % save power envelopes
        if Settings.SaveCorrected.envelopes,
            saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
            ROInets.make_directory(saveDir);
            saveFile = fullfile(saveDir,                                      ...
                                sprintf('%s_%s_ROI_envelope_timecourses.mat', ...
                                        sessionName, bandName));
            save(saveFile, 'nodeEnv');
        end%if
        
        % calculate correlation matrices. 
        CorrMats{iFreq} = ROInets.run_correlation_analysis(nodeData,     ...
                                                      nodeEnv,      ...
                                                      Settings.Regularize);
        % free up some memory
        clear nodeData nodeEnv
    end%if
    
    %% Use an empirical null to enable conversion to z-stats
    sigma = ROInets.find_empirical_H0_distribution_width(nH0Iter,                   ...
                                                         nNodes,                    ...
                                                         nTimeSamples,              ...
                                                         Settings,                  ...
                                                         CorrMats{iFreq}.Regularization, ...
                                                         ARmodel,                   ...
                                                         Fs, Filter,                ...
                                                         Settings.EnvelopeParams);
                                                     
    %% store AR model parameters
    CorrMats{iFreq}.ARmodel = ARmodel;
    CorrMats{iFreq}.H0Sigma = sigma;
    
    %% Store session name
    CorrMats{iFreq}.sessionName = sessionName;
    CorrMats{iFreq}.timeWindow  = timeRange;
    
    %% conversion of correlations to z-stats
    fprintf(' Converting correlations to normal z-stats\n');
    CorrMats{iFreq} = ROInets.convert_correlations_to_normal_variables(CorrMats{iFreq}, ...
                                                                  sigma,      ...
                                                                  doRegularize);
    
    % average correlations from pairwise analysis - best to average after
    % conversion to z-stats. 
    if strcmpi(Settings.leakageCorrectionMethod, 'pairwise'),
        CorrMats{iFreq} = ROInets.apply_function_to_structure(@balance_correlations, ...
                                                         CorrMats{iFreq});
    end%if
end%loop over freq bands

%% save results to disc to provide backup of completed session analyses
save(resultsSaveName, 'CorrMats');




%%% END OF FUNCTION PROPER %%%














%% SUBFUNCTIONS
% use a nested function to minimise memory movement
%--------------------------------------------------------------------------
function [Filter, bandName] = filter_voxel_data(Settings, Fs, iFreq)
%FILTER_VOXEL_DATA band-pass filters voxel data
%
% [FILTER, BANDNAME] = FILTER_VOXEL_DATA(SETTINGS, FS, IFREQ)
%   returns the filtered data VOXELDATA, the FILTER settings used and
%   BANDNAME, a string with the frequency band name. SETTINGS is a
%   structure returned from ROINETS.CHECK_INPUTS, Fs is the sampling
%   rate and iFREQ is the index of the frequency bands to use. 
%
% The function is nested, running on variable VOXELDATABROADBAND and
% returning VOXELDATA. 

if isempty(Settings.frequencyBands{iFreq}),
    fprintf(' No filtering applied\n');

    bandName  = 'broadBand';
    voxelData = voxelDataBroadBand;
    Filter    = struct('band', []);
else
    fprintf(' Filtering for band %f-%f Hz\n',  ...
            Settings.frequencyBands{iFreq}(1), ...
            Settings.frequencyBands{iFreq}(2));

    bandName  = sprintf('%d-%dHz', ...
                       Settings.frequencyBands{iFreq}(1), ...
                       Settings.frequencyBands{iFreq}(2));
    Filter    = struct('order',     4,         ...
                       'type',      'but',     ...  % butterworth IIR filter
                       'direction', 'twopass', ...  % no phase lag or edge effects
                       'band',      Settings.frequencyBands{iFreq});

    voxelData = ft_preproc_bandpassfilter(voxelDataBroadBand, ...
                                          Fs,                 ...
                                          Filter.band,        ...
                                          Filter.order,       ...
                                          Filter.type,        ...
                                          Filter.direction);
end%if
end%filter_voxel_data










%--------------------------------------------------------------------------
function [nodeData, voxelWeightings] = get_corrected_node_tcs(spatialBasis, ...
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
%     'PCA'   - take 1st PC of variance-normalised voxels 
%     'mean'  - mean voxel time-course.
%     'peakVoxel' - voxel with the maximum variance
%
% NODEDATA = GET_CORRECTED_NODE_TCS(VOXELDATA, SPATIALBASIS, 'spatialBasis', METHOD)
%   it is also possible to use a spatial basis set (e.g. from group ICA) to
%   infer parcel time-courses. Each spatial map (held in columns) is a
%   whole-brain map - each map can be non-binary and maps may overlap.
%   SPATIALBASIS is therefore an (nVoxels x nSamples) matrix. The ROI
%   time-course for each spatial map is the 1st PC from all
%   variance-normalised voxels, weighted by the spatial map. 
%
% [NODEDATA, VOXEL_WEIGHTINGS] = GET_CORRECTED_NODE_TCS(VOXELDATA, SPATIALBASIS, 'spatialBasis', METHOD)
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
%	$Revision: 263 $
%	$LastChangedDate: 2014-10-23 11:30:39 +0100 (Thu, 23 Oct 2014) $
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
                              
                voxelWeightings(thisMask, iParcel) = U ./ sum(abs(U));
                
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
        
        
    case 'mean'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'It will be binarised. \n']);
        end%if
        spatialBasis = logical(spatialBasis);
        
        % take mean time-course in each parcel
        for iParcel = nParcels:-1:1,
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding mean time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            thisMask = spatialBasis(:, iParcel);
            if any(thisMask), % non-zero
                parcelData              = voxelData(thisMask, :);
                nodeDataOrig(iParcel,:) = nanmean(parcelData, 1);          % mean at each point in time. Produces a row vector.
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
        if nargout > 1, 
            voxelWeightings = spatialBasis;
        end%if
        
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
        
        % find time-course for each spatial basis map
        for iParcel = nParcels:-1:1, % allocate memory on the fly
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding spatial basis time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            % extract the spatial map of interest
            thisMap     = scaledSpatialMaps(:, iParcel);
            
            % estimate temporal-STD for normalisation
            temporalSTD = max(std(voxelData, [], 2), eps);
            
            % variance-normalise all voxels to remove influence of
            % outliers. - remove this step 20 May 2014 for MEG
            % Weight all voxels by the spatial map in question
            weightedTS  = ROInets.scale_rows(voxelData, thisMap);
            
            % perform svd and take scores of 1st PC as the node time-series
            % U is nVoxels by nComponents - the basis transformation
            % S*V holds nComponents by time sets of PCA scores - the 
            % timeseries data in the new basis
            [U, S, V]   = ROInets.fast_svds(weightedTS, 1);
            clear weightedTS
            
            PCAscores   = S * V';
            maskThresh  = 0.5; % 0.5 is a decent arbitrary threshold chosen by Steve Smith and MJ after playing with various maps.
            thisMask    = thisMap > maskThresh;   
            
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
                      
                voxelWeightings(:, iParcel)         = U ./ ...
                                                      sum(abs(U(thisMask)));
                voxelWeightings(~thisMask, iParcel) = 0;
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(PCAscores));
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
%--------------------------------------------------------------------------
end%run_individual_correlation_analysis
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
function balanced = balance_correlations(x)
%BALANCE_CORRELATIONS makes matrices symmetric

if ismatrix(x) && ~isvector(x),
    balanced = (x + x') ./ 2.0;
else
    balanced = x;
end%if
end%balance_correlations






%--------------------------------------------------------------------------
function [t, tI, tR] = time_range(time, timeRange, iSession)
%TIME_RANGE selects time range for analysis of each session
% TIME is a vector of times
% TIMERANGE is either a cell array of two-component vectors, a single
% two-component vector, or a null vector

if isempty(timeRange),
    % use the whole time range
    t = time;
    tI = true(size(time));
    tR = [];
else
    % subselect time range
    if iscell(timeRange),
        tR = timeRange{iSession};
    else
        tR = timeRange;
    end%if
    validateattributes(tR, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nonnegative', 'nondecreasing'}, ...
                       'time_range', 'timeRange', 2);
                   
    tI = (time <= tR(2)) & (time >= tR(1));
    t  = time(tI);
end%if
end%time_range





%--------------------------------------------------------------------------
function [] = save_corrected_timecourse_results(nodeData, ...
    allROImask, voxelWeightings, Settings, sessionName, protocol, bandName)
% Save the weightings over voxels used to calculate the time-course
% for each ROI
if Settings.SaveCorrected.ROIweightings,
    allVoxelWeightings                = zeros(ROInets.rows(allROImask), ...
                                              ROInets.cols(voxelWeightings));
    allVoxelWeightings(allROImask, :) = voxelWeightings;

    saveDir             = fullfile(Settings.outputDirectory, ...
                                   'spatialBasis-ROI-weightings', filesep);
    ROItcWeightSaveFile = fullfile(saveDir,                      ...
                                   [sessionName '_' protocol '_' ...
                                    bandName '_ROI_timecourse_weightings']);
    ROInets.make_directory(saveDir);
    nii_quicksave(allVoxelWeightings, ROItcWeightSaveFile, Settings.gridStep);
end%if

% save node data
if Settings.SaveCorrected.timeCourses,
    saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
    ROInets.make_directory(saveDir);
    saveFile = fullfile(saveDir, [sessionName '_correction-' protocol '_' bandName '_ROI_timecourses.mat']);
    save(saveFile, 'nodeData');
end%if

% save variance in each ROI
if Settings.SaveCorrected.variances,
    saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
    ROInets.make_directory(saveDir);
    varSaveFile = fullfile(saveDir, [sessionName '_correction-' protocol '_' bandName '_ROI_variances.mat']);
    ROIvariances = var(nodeData, [], 2);                                   %#ok<NASGU>
    save(varSaveFile, 'ROIvariances');
end%if
end%save_corrected_timecourse_results
%--------------------------------------------------------------------------
% [EOF]