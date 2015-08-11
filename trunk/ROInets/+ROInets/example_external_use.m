function correlationMats = example_external_use(varargin)
%EXAMPLE  example analysis for a single subject, if you're external to OHBA
%
% correlationMats = EXAMPLE_EXTERNAL_USE()
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
%	Originally written on: MACI64 by Giles Colclough, 11-Aug-2015 


% set path to source-reconstructed data
dataFile = '/path/to/data/source_data_file.mat';

% create an toy dataset
% use an autoregressive model to generate smooth data
% follow http://www.mathworks.co.uk/help/signal/examples/linear-prediction-and-autoregressive-modeling.html
Fs       = 100; %Hz
duration = 60; %s
time     = 0:1.0/Fs:duration;
nSamples = length(time);
b        = fir1(1024, 0.5);
nVoxels  = 3;
[ARfilterTerms, ARnoiseVar] = lpc(b, 7); 
% generate data from a covariance matrix and smooth
C    = [1  -0.1 0.6
       -0.1   1 0.3
        0.6 0.3   1] * ARnoiseVar;
u    = chol(C)' * randn(nVoxels, nSamples);
data = filter(1, ARfilterTerms, u.').';
figure('Name', 'Input data', 'Color', 'w');
plot(time.', data.');
% save to file
sampleRateInHz = Fs;
save(dataFile, 'data', 'time', 'sampleRateInHz');

% choose a binary ROI map. 
spatialBasis = eye(3);
                  
% choose a results directory
outDir = '/path/to/results/';

% set a save file name
resultsName = fullfile(outDir, 'myResults');
sessionName = 'myBestSubject';

% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = spatialBasis;                   % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
Settings.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
Settings.frequencyBands           = {[13 30]};                      % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = sessionName; 
Settings.SaveCorrected            = struct('timeCourses',   false, ...  % save corrected timecourses
                                           'envelopes',     true,  ...  % save corrected power envelopes
                                           'variances',     false);     % save mean power in each ROI before correction

% run the ROI network analysis
Settings        = ROInets.check_inputs(Settings);
correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);

% Want to run an analysis on many subjects? Have a look at 
% osl_network_analysis to see the suggested steps. 

% show results
figure('Name', 'node correlation matrix', 'Color', 'w');
imagesc(correlationMats{1}.envCorrelation);
axis square
colorbar

% plot envelopes
load(fullfile(outDir, 'corrected-ROI-timecourses', ...
              'myBestSubject_13-30Hz_ROI_envelope_timecourses.mat'));
figure('Name', 'envelope timecourses', 'Color', 'w');
plot(time_ds', nodeEnv')
end%example_1subj
% [EOF]