function correlationMats = osl_network_analysis(Dlist, Settings)
%OSL_NETWORK_ANALYSIS	ROI correlation and partial correlation matrices
%
% [CORRELATIONMATS] = osl_network_analysis(DLIST, Settings) runs a
%   correlation-based network analysis between ROIs on a set of
%   source-reconstructed data. These are provided as a cell array of SPM12
%   MEEG objects, or their filenames, in DLIST. 
%
%   Relevant input parameters are held in Settings:
%      spatialBasisSet          - the name of a binary nifti file which 
%                                 holds the voxel allocation for each ROI
%                                 (each ROI allocation is a volume), or a
%                                 spatial basis set from ICA (each spatial
%                                 map is a volume)
%
%      subjectsToDo             - vector of indices indicating which of the 
%                                 source reconstruction sessions to analyse
%
%      nParallelWorkers         - analyse several sessions in parallel, 
%                                 if non-zero. 
%                                 Specifies number of matlab workers to use. 
%                                 Caution: these tasks are very memory-heavy! 
%                                 You might slow down because you max out 
%                                 on RAM. 
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
%      leakageCorrectionMethod  - choose from 'symmetric', 'pairwise' or 
%                                 'none'.
%                                 The symmetric method is recommended. It
%                                 orthogonalises the ROI time-courses in an
%                                 all-to-all symmetric fashion. 
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
%      envelopeWindowLength     - sliding window length for power envelope 
%                                 calculation. 2 s is a good value. 
%                                 See Brookes 2011, 2012 and Luckhoo 2012. 
%
%      envelopeWindowOverlap    - Overlap on the sliding window. Try 0.6?   
%
%      envelopeWithFilter       - Set true or false to use a more 
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
%      timecourseCreationMethod - 'PCA' or 'peakVoxel'  for binary 
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
%     fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses')
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
%       envCorrelation_z                        : z-stats for envelope
%                                                 correlations
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelation_z                 : z-stats for envelope
%                                                 partial correlations
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelationRegularized_z      : z-stats for envelope
%                                                 partial correlaitons
%                                                 nROIs x nROIs x nSessions
%
%       groupEnvCorrelation_z                   : group result of z-test on
%                                                 mean of first-level
%                                                 z-statistics
%
%       groupEnvPartialCorrelation_z            : group result of z-test on
%                                                 mean of first-level
%                                                 z-statistics
%
%       groupEnvPartialCorrelationRegularized_z : group result of z-test on
%                                                 mean of first-level
%                                                 z-statistics
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
%   References:
%   Colclough, G.L. et al., "A symmetric multivariate leakage correction
%   for MEG connectomes," NeuroImage 117, pp. 439-448 (2015). 



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


%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 17:53:55

fprintf('%s: starting analysis. \n', mfilename);

%% Parse input settings
fprintf('%s: checking inputs. \n', mfilename);
Settings = ROInets.check_inputs(Settings);

assert(numel(Dlist) == Settings.nSessions, ...
       [mfilename ':WrongNoSessions'],     ...
       'Number of sessions should match number of D objects passed in. \n');

% make results directory
ROInets.make_directory(Settings.outputDirectory);

% save settings
outputDirectory = Settings.outputDirectory;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

%% Run correlation analysis on each subject

fprintf('%s: Running correlation analysis. \n', mfilename);

for iSession = Settings.nSessions:-1:1,
    fprintf('\n\n%s: Individual correlation analysis for file %d out of %d\n', ...
            mfilename, Settings.nSessions - iSession + 1, Settings.nSessions);

    D                          = Dlist{iSession};
    sessionName                = Settings.sessionName{iSession};
    matsSaveFileName{iSession} = fullfile(outputDirectory,                                      ...
                                          sprintf('%s_single_session_correlation_mats_tmp.mat', ...
                                                  sessionName));

    mats{iSession} = ROInets.run_individual_network_analysis(D,                          ...
                                                             Settings,                   ...
                                                             matsSaveFileName{iSession}, ...
                                                             iSession);
end%for

% reformat results - correlationMats is a cell array of frequency bands
correlationMats = ROInets.reformat_results(mats, Settings);
clear mats

%% Group-level analysis
correlationMats = ROInets.do_group_level_statistics(correlationMats, Settings);

%% save matrices
fprintf('\n%s: Saving Results. \n', mfilename);

% save collected results
saveFileName = fullfile(outputDirectory, 'ROInetworks_correlation_mats.mat');
save(saveFileName, 'correlationMats');

% we stored individual results as we went along, in case of crash. Delete
% them if we've safely got to this stage. 
for iSession = 1:length(matsSaveFileName),
    delete(matsSaveFileName{iSession});
end%for
 
% tidy output of funciton
Settings.correlationMatsFile = saveFileName;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

fprintf('%s: Analysis complete. \n\n\n', mfilename);
end%osl_network_analysis
% [EOF]