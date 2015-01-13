function FAIL = test_me()
%TEST_ME  test that pipeline runs

%	Copyright 2015 OHBA
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


%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@magd.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 12-Jan-2015 15:17:08

% These are hardly comprehensive unit tests, but better than nowt...


%% General things
tilde = '/Users/gilesc';

% set path to source-reconstructed data
dataDir = fullfile(tilde, 'data', 'ROInets-test', filesep);
ROInets.make_directory(dataDir);

tmp = load('/Users/gilesc/data/Henry/OIL-output/analysis-for-parcellation-paper/OIL_4-30Hz_8mm.oil/oil_2s_win_8mm_grid_concatenated_subject_data_spatial_ica_40comps_0_its/oil.mat', 'oil');
[data, time, sampleRateInHz] = MEGsim.get_source_data(tmp.oil, '9_ctf_3-state_eo_recon');
dataFile = fullfile(dataDir, 'source_data.mat');
save(dataFile, 'data', 'time', 'sampleRateInHz', '-v7.3');
clear data time Fs

% choose a binary ROI map. Take care that the resolution of the nifti file
% matches that of the source reconstruction.
fMRIparcellationDirectory = '/Users/gilesc/data/papers/parcellation-orthogonalisation/parcellations/fMRI-parcellations/groupPCA_d1000_d100.dr';
parcelFile = fullfile(fMRIparcellationDirectory, ...
                      'fmri_d100_parcellation_reduced_2mm_ds8mm.nii.gz');
                  
% choose a results directory
outDir = dataDir;

% set a save file name
resultsName = fullfile(outDir, 'testResults');
sessionName = 'myBestSubject';


%% TEST 1
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = {[60 120]};                             % range of times to use for analysis
Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
Settings.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
Settings.frequencyBands           = {[8 13]};                       % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelope      = true;
Settings.SaveCorrected.variances     = true;
Settings.SaveCorrected.ROIweightings = true;

% run the ROI network analysis
try
    correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);
catch ME
    fprintf('%s: Test 1 failed. \n', mfilename);
    rethrow(ME);
end%try

%% TEST 2
binaryparcellationDirectory = '/Users/gilesc/data/papers/parcellation-orthogonalisation/parcellations/parcellation-for-sim';
parcelFile = fullfile(binaryparcellationDirectory, ...
                      'parcellation_for_simulation_8_mm_Aug.nii.gz');
 
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = [];                             % range of times to use for analysis
Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.Regularize.method        = 'Bayesian';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
Settings.Regularize.Prior         = struct('a', 1./3, 'b', 0);
Settings.leakageCorrectionMethod  = 'symmetric';                      % choose from 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
Settings.frequencyBands           = {[8 13]};                       % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'PCA';                          % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelope      = false;
Settings.SaveCorrected.variances     = false;
Settings.SaveCorrected.ROIweightings = false;

% run the ROI network analysis
try
    correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);
catch ME
    fprintf('%s: Test 2 failed. \n', mfilename);
    rethrow(ME);
end%try

%% TEST 3
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = {[0 120]};                             % range of times to use for analysis
Settings.Regularize.do            = false;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.leakageCorrectionMethod  = 'pairwise';                      % choose from 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.overlap   = 0.5;
Settings.EnvelopeParams.useFilter    = false;                        % use a more sophisticated filter than a sliding window average
Settings.frequencyBands           = {[]};                       % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'mean';                          % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelope      = false;
Settings.SaveCorrected.variances     = false;
Settings.SaveCorrected.ROIweightings = false;

% run the ROI network analysis
try
    correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);
catch ME
    fprintf('%s: Test 3 failed. \n', mfilename);
    rethrow(ME);
end%try

%% We've got to the end!
FAIL = 0;
ROInets.call_fsl_wrapper(['rm -rf ' dataDir]);

end%test_me
% [EOF]