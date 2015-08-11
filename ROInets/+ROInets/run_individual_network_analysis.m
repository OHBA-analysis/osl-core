function CorrMats = run_individual_network_analysis(D,        ...
                                                    Settings, ...
                                                    resultsSaveName, ...
                                                    iSession)
%RUN_INDIVIDUAL_CORRELATION_ANALYSIS runs a single session network analysis
%
% CORRMATS = RUN_INDIVIDUAL_NETWORK_ANALYSIS(SPMMEEG, SETTINGS, SAVENAME)
%   performs a network correlation analysis between ROIs on a single subject. 
%
%   SPMMEEG is an SPM MEEG object.
%     data: nvoxels x time MEG data for analysis for a single subject. It
%           is assumed that the number of voxels matches the number of voxels in
%           the ROI basis set, which is passed in as a nifti file (see below). 
%
%
%   RESULTS_SAVE_NAME sets the location where the individual session
%   correlation matrices will be saved to disc. 
%
%
% CORRMATS = RUN_INDIVIDUAL_NETWORK_ANALYSIS(VOXELDATAFILE, SETTINGS, SAVENAME)
%   uses the name of a .mat file on disk, VOXELDATAILE. This .mat file
%   should contain three variables:
%     data: nvoxels x time MEG data for analysis for a single subject. It
%           is assumed that the number of voxels matches the number of voxels in
%           the ROI basis set. 
%
%     time: time-course for the MEG data
%     sampleRateInHz: also required. A constant sampling rate is assumed.
%
%   Relevant input parameters are held in SETTINGS:
%      spatialBasisSet          - the name of a binary nifti file which 
%                                 holds the voxel allocation for each ROI
%                                 (each ROI allocation is a volume), or a
%                                 spatial basis set from ICA (each spatial
%                                 map is a volume)
%                               - alternatively, pass in a voxels x ROIs
%                                 matrix. 
%
%      gridStep                 - grid step in mm of spatialBasisSet nifti
%                                 file. 
%
%      timeRange                - two-component vector of times to
%                                 sub-select a range for analysis. Leave as
%                                 empty ([]) to use the entire range
%                                 (default behaviour). Also can be passed
%                                 as a cell array of two-component vectors,
%                                 one for each subject. 
%      
%      Regularize               - Structure controlling the use of
%                                 regularization to find partial
%                                 correlation matrices, using the graphical
%                                 lasso. It has fields:
%            .do            : true or false: controls use of
%                             regularization. [false]
%            .method        : 'Bayesian' or 'Friedman': use Wang (2012)'s
%                             Bayesian graphical lasso or Friedman (2007)'s
%                             graphical lasso
%            .path          : This specifies a single, or vector of, 
%                             possible rho-parameters controlling the
%                             strength of regularization.
%                             If a vector is selected, 10-fold cross-validation
%                             is used to pick the optimal rho-parameter
%                             for each subject. The corrected Akaike
%                             Information Criterion is used to select
%                             the optimal rho. 
%                             [Friedman only]
%            .adaptivePath  : true or false: adapt path if the best
%                             regularization parameter is on the edge of
%                             the path. Then finesse the path three times. 
%                             [Friedman only] 
%            .Prior         : structure with fields controlling the shape
%                             of the gamma(x; a, 1/b) hyperprior on the 
%                             regularization parameter. 
%                             If this field is not set, the default Kerman 
%                             (2011) neutral hyperprior is used (a=1/3, b=0)
%                             [Bayesian only]
%                .a - shape
%                .b - 1/scale
%           for references, type `help ROInets.run_correlation_analysis'. 
%
%           If the Bayesian method is selected and the nEmpiricalSamples
%           options is used (see below), each empirical sample will use a
%           different regularization parameter at random from the posterior
%           distribution of regularization parameters used as part of the
%           MCMC scheme. nEmpiricalSamples should be relatively high to
%           adequately sample the space of 8000 values. 
%
%      leakageCorrectionMethod  - choose from 'symmetric', 'closest', 
%                                 'pairwise' or 'none'.
%                                 The symmetric method orthogonormalises 
%                                 the ROI time-courses in an
%                                 all-to-all symmetric fashion. 
%                                 The closest method starts with the
%                                 symmetric method, then lifts the
%                                 normality constraint to iterate to the
%                                 closest orthogonal matrix. 
%                                 The pairwise method extends voxel-wise
%                                 leakage correction methods to the parcel
%                                 time-courses. Correlations between pairs
%                                 of nodes are the average of correlations
%                                 found when one node is regressed from the
%                                 other, and then visa-versa.
%                                 'None' applies no spatial leakage
%                                 correction. 
%
%      nEmpiricalSamples        - convert correlations to standard normal 
%                                 z-statistics using a simulated empirical 
%                                 distribution. This controls how many 
%                                 times we simulate null data of the same 
%                                 size as the analysis dataset. 
%                                 Set to zero to make normal assumptions, 
%                                 and decrease runtime. 
%                                 For datasets with many ROIs, only a few
%                                 samples are required - O(10). 
%
%      ARmodelOrder             - We tailor the empirical data to have the 
%                                 same temporal smoothness as the MEG data.
%                                 Set an order of 0 to use normality
%                                 assumptions. 
%                                 An order of 1 should be ok, but you could
%                                 use partial autocorrelation on some of 
%                                 the voxel time-courses to check this 
%                                 assumption. (Check out the matlab function
%                                 aryule, and google for partial
%                                 autocorrelation.)
%
%      EnvelopeParams           - Structure controlling how time-courses
%                                 are enveloped and down-sampled. Has
%                                 fields:
%
%          .windowLength :        sliding window length for power envelope 
%                                 calculation. 2 s is a good value. 
%                                 See Brookes 2011, 2012 and Luckhoo 2012. 
%
%          .overlap      :        Overlap on the sliding window. Try 0.6?  
%
%          .useFilter    :        Set true or false to use a more 
%                                 sophisticated downsamping operation than 
%                                 a moving average to find the power envelope. 
%                                 The passed frequency is 1.0/envelopeWindowLength.
%                                 The window overlap parameter is then
%                                 irrelevant. 
%
%      frequencyBands           - cell array of frequency bands for analysis. 
%                                 Set to empty to use broadband. 
%                                 E.g. {[4 8], [8 13], []}
%                                 The bandpass filtering is performed 
%                                 before orthogonalisation. 
%
%      timecourseCreationMethod - 'PCA', 'peakVoxel' or 'mean' for binary 
%                                 parcellations. 
%                                 Sets whether an ROI time-course is 
%                                 created using the mean over all voxels in 
%                                 the ROI, or by choosing the coefficients 
%                                 of the principal component accounting for 
%                                 the majority of the ROI's variance, or by
%                                 selecting the voxel with the greatest
%                                 variance.
%
%                                 If a spatial basis set (e.g. from ICA) is
%                                 passed in, 'spatialBasis' can be set,
%                                 which uses the PCA method, and
%                                 incorporating the whole-brain weightings
%                                 of each spatial map. 
%
%      groupStatisticsMethod    - Choose 'mixed-effects' or 'fixed-effects'
%                                 for group-level analysis. The former
%                                 performs a t-test on the means of the
%                                 first-level z-stats; the latter performs
%                                 a z-test. The t-test may have less power
%                                 but is also less sensitive to
%                                 inaccuracies in the scaling of
%                                 first-level z-scores - it is preferred,
%                                 but may need a large number of subjects
%                                 to produce interpretable results.
%
%      outputDirectory          - Set a directory for the results output
%
%      sessionName              - string used to construct filenames for
%                                 saving outputs. 
%
%      SaveCorrected            - Structure determining what is saved to
%                                 disc. Has fields:
%
%          .timeCourses     : orthgonalised ROI time-courses
%          .envelopes       : corrected ROI envelopes
%          .variances       : variances in each ROI
%          .ROIweightings   : weights used in each ROI used to construct
%                             single time-course
%
%
%   The ROI time-courses, corrected for source spread, and their envelopes, 
%   are saved in 
%     fullfile(oil.ROInetworks.outputDirectory, 'corrected-ROI-timecourses')
%
%   Various correlation matrices are output in CORRELATIONMATS. The
%   variable is a cell array, with one cell for each frequency band in the
%   analysis. Each cell holds a structure with fields:
%      
%       correlation                             : correlation between 
%                                                 corrected time-courses 
%                                                 (expected to be null) 
%                                                 nROIs x nROIs x nSessions
%
%       envCorrelation                          : correlation between
%                                                 corrected BLP envelopes
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelation                   : partial correlation
%                                                 between BLP envelopes
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelationRegularized        : L1-regularized partial 
%                                                 correlation between BLP 
%                                                 envelopes
%
%       env_z                                   : z-stats for envelope
%                                                 correlations
%                                                 nROIs x nROIs x nSessions
%
%       env_z_partial                           : z-stats for envelope
%                                                 partial correlations
%                                                 nROIs x nROIs x nSessions
%
%       env_z_partial_reg                       : z-stats for envelope
%                                                 partial correlations
%                                                 nROIs x nROIs x nSessions
%
%       Regularization                          : Chosen rho-parameter for
%                                                 L1 regularization in each
%                                                 subject
%
%       ARmodel                                 : AR coefficients estimated
%                                                 from MEG data for each
%                                                 session
%
%       H0sigma                                 : Width of empirical H0
%                                                 null correlation z-stats
%
%   These are also saved to disk at: 
%   fullfile(Settings.outputDirectory, 'ROInetworks_correlation_mats.mat')
%   
%   The analysis settings are also saved in the output directory. 
%
%   The methodology used in this pipeline is set out in 
%   Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A
%   symmetric multivariate leakage correction for MEG connectomes,"
%   NeuroImage [in revision]. 



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

fprintf('%s: starting analysis. \n', mfilename);


%% Parse input settings
fprintf('%s: checking inputs. \n', mfilename);

% Settings = ROInets.check_inputs(Settings);

% make results directory
ROInets.make_directory(Settings.outputDirectory);

% save settings
outputDirectory = Settings.outputDirectory;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

% load parcellation
parcelFile = Settings.spatialBasisSet;
if ~isempty(parcelFile) && ischar(parcelFile), 
    spatialBasis = nii_quickread(parcelFile, Settings.gridStep);
elseif (islogical(parcelFile) || isnumeric(parcelFile)) && ismatrix(parcelFile), % expecting vox x parcels
    spatialBasis = parcelFile;
else
    error([mfilename ':ParcelReadFail'], ...
          'Failed to read parcel file %s. \n', ...
          parcelFile);
end%if
clear parcelFile

% is this function being run as part of a loop?
if nargin < 4 || ~exist('iSession', 'var') || isempty(iSession),
    iSession    = 1;
    sessionName = Settings.sessionName;
else
    sessionName = Settings.sessionName{iSession};
end%if
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
% load in data
try
    % assume spm meeg object
    % parse any string inputs automatically
    D                  = spm_eeg_load(D); 
    fprintf('%s: loading data from file %s. \n', mfilename, D.fname);
    allTime            = D.time;
    Fs                 = D.fsample;
    badSamples         = any(D.badsamples(:,:,:));
    voxelDataBroadBand = D(:,find(~badSamples),:);                         %#ok<FNDSB>

catch ME
    try % what about the name of a .mat file?
        assert(2 == exist(D, 'file'),                           ...
               [mfilename ':UnrecognisedDataInput'],            ...
               ['The first input should be an spm MEEG object', ...
                'or path to a .mat file. \n']);
        fprintf('%s: loading data from file %s. \n', mfilename, D);
        tmpInput           = load(D);
        voxelDataBroadBand = tmpInput.data;
        Fs                 = tmpInput.sampleRateInHz;                      
        allTime            = tmpInput.time;
        badSamples         = false(size(allTime));                         
        clear tmpInput
    catch
        rethrow(ME);
    end%try
end%try


% select time range from user-defined input
[time, timeInd, timeRange] = time_range(allTime(~badSamples), ...
                                        Settings.timeRange, iSession);

% apply to source data, too.
voxelDataBroadBand(:,~timeInd) = [];
nTimeSamples                   = ROInets.cols(voxelDataBroadBand);
nNodes                         = ROInets.cols(spatialBasis);

assert(sum(timeInd) == nTimeSamples,           ...
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
    

    [nodeDataUnCorr, voxelWeightings] = get_node_tcs(real(spatialBasis), ...
                                                     tcsMethod);
    nodeData = ROInets.remove_source_leakage(nodeDataUnCorr, protocol);
    save_corrected_timecourse_results(nodeData, allROImask,               ...
                                      voxelWeightings, Settings,          ...
                                      sessionName, protocol, bandName);
    clear nodeDataUnCorr;
    
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
        [nodeEnv, time_ds] = ROInets.envelope_data(nodeData,        ...    
                                                   time,            ...
                                                   Settings.EnvelopeParams); %#ok<ASGLU>
        
        % save power envelopes
        if Settings.SaveCorrected.envelopes,
            saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
            ROInets.make_directory(saveDir);
            saveFile = fullfile(saveDir,                                      ...
                                sprintf('%s_%s_ROI_envelope_timecourses.mat', ...
                                        sessionName, bandName));
            save(saveFile, 'nodeEnv', 'time_ds');
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
function [nodeData, voxelWeightings] = get_node_tcs(spatialBasis, ...
                                                    timeCourseGenMethod)
%GET_NODE_TCS correct ROI time-courses for source leakage
%
% methods for extracting node time courses using differing
% orthogonalisation protocols
%
% NODEDATA = GET_NODE_TCS(VOXELDATA, PARCELFLAG, METHOD)
%   produces orthogonalised node time-courses NODEDATA from 
%   (nVoxels x nSamples) matrix of beamformed data VOXELDATA. VOXELDATA may
%   be passed in as an array, or as the filename of a saved .MAT file. 
%
%   PARCELFLAG is a logical array (nVoxels x nParcels) which identifies
%   the voxels in each parcel. Each column identifies the voxels 
%   making up each parcel. True entries indicate membership. 
%
%   METHOD chooses how the ROI time-course is created. It can be:
%     'PCA'   - take 1st PC of voxels 
%     'peakVoxel' - voxel with the maximum variance
%     note that the mean timecourse suffers issues relating to sign
%     ambiguities, and so is not currently an option.
%
% NODEDATA = GET_NODE_TCS(VOXELDATA, SPATIALBASIS, METHOD)
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
% [NODEDATA, VOXEL_WEIGHTINGS] = GET_NODE_TCS(...)
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
            
            nodeData(iParcel,:) = nodeTS;
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
                nodeData(iParcel,:) = voxelData(maxPowerInd,:);
                
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
                     
                nodeData(iParcel,:) = zeros(1, ROInets.cols(voxelData));
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
        
        % find time-course for each spatial basis map
        for iParcel = nParcels:-1:1, % allocate memory on the fly
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [' ' mfilename ...
                         ':    Finding spatial basis time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            % extract the spatial map of interest
            thisMap     = scaledSpatialMaps(:, iParcel);
            parcelMask  = logical(thisMap);
            
            % estimate temporal-STD for normalisation
            temporalSTD = max(std(voxelData, [], 2), eps);
            
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
            
            nodeData(iParcel, :) = nodeTS;
            
        end%loop over parcels
        
        clear voxelData 
        
        
    otherwise
        error([mfilename ':UnrecognisedTimeCourseMethod'],            ...
              ['Unrecognised method for finding ROI time-course. \n', ...
               'Expected ''PCA'', ''spatialBasis'', or ''mean''. \n']);
end%switch

ft_progress('close');

end%get_node_tcs
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
    try
        nii_quicksave(allVoxelWeightings, ROItcWeightSaveFile, Settings.gridStep);
    catch % perhaps we have a weird number of voxels
        save(ROItcWeightSaveFile, 'allVoxelWeightings');
    end%try
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