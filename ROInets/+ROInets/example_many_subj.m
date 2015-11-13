function correlationMats = example_many_subj(varargin)
%EXAMPLE  example analysis for a single subject
%
% correlationMats = EXAMPLE_MANY_SUBJ()
% for help, type
% `help ROInets.run_individual_network_analysis'

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
%	Originally written on: MACI64 by Giles Colclough, 12-Jan-2015 15:25:47


% set path to source-reconstructed data
dataDir = '/path/to/data/';

% choose a binary ROI map. Take care that the resolution of the nifti file
% matches that of the source reconstruction.
parcellationDirectory = '/path/to/parcellation/';
parcelFile = fullfile(parcellationDirectory, ...
                      'parcellationFile.nii.gz');
                  
% choose a results directory
outDir = '/path/to/results/';

% set objects in and names for each session
for iFile = 1:10; % loop over D objects in directory
    Dlist{iFile}       = fullfile(dataDir, ...
                                  sprintf('source_data_SPM_obj_%d.mat', iFile));
    sessionName{iFile} = sprintf('obj_%d', iFile);
end%for

% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
Settings.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
Settings.frequencyBands           = {[8 13], [13 30], []};          % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected            = struct('timeCourses',   false, ...  % save corrected timecourses
                                           'envelopes',     true,  ...  % save corrected power envelopes
                                           'variances',     false);     % save mean power in each ROI before correction

% run the ROI network analysis
correlationMats = osl_network_analysis(Dlist, Settings);
end%example_many_subj
% [EOF]