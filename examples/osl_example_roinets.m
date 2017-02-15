%% Example script for ROInets-based analysis

osldir = getenv('OSLDIR');
data_dir = fullfile(osldir,'example_data','roinets_example');
output_directory = fullfile(osldir,'practical','roinets_demo');
mkdir(output_directory)

%%
% Copy the demo files and make a list of the D files
subjects = 1:10;
D_files = {};
session_name = {};
mkdir(fullfile(osldir,'practical','roinets_demo'));
for j = 1:length(subjects)
    session_name{j} = sprintf('subject_%d',j);
    D = spm_eeg_load(fullfile(data_dir,session_name{j}));
    D_files{j} = D.copy(fullfile(output_directory,session_name{j}));
    D_files{j} = D_files{j}.montage('switch',2);
end

%% If files already exist
subjects = 1:10;
D_files = {};
session_name = {};
for j = 1:length(subjects)
    session_name{j} = sprintf('subject_%d',j);
    D_files{j} = spm_eeg_load(fullfile(data_dir,session_name{j}));
    D_files{j} = D_files{j}.montage('switch',2);
end


%% 
% Load the parcellation and save a parcelflag variable
% Settings.spatialBasisSet = fullfile(output_directory,'parcelflag.nii.gz');
% p = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
% p.savenii(p.to_matrix(p.binarize),Settings.spatialBasisSet);
% Settings.gridStep = p.resolution;

%%
% setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'); % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
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
Settings.EnvelopeParams.takeLogs  = true;                           % perform analysis on logarithm of envelope. This improves normality assumption
Settings.frequencyBands           = {[8 13], [13 30], []};          % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA',  'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = output_directory;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = session_name; 
Settings.SaveCorrected            = struct('timeCourses',   false, ...  % save corrected timecourses
                                           'envelopes',     true,  ...  % save corrected power envelopes
                                           'variances',     false);     % save mean power in each ROI before correction

%%
% Could call osl_network_analysis at this point
ROInets.osl_network_analysis(D_files,Settings)

%% 
% Or do the same manuall



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

    if strcmpi(Settings.paradigm, 'task'),
        mats{iSession} = ROInets.run_individual_network_analysis_task(D,                          ...
                                                                 Settings,                   ...
                                                                 matsSaveFileName{iSession}, ...
                                                                 iSession);
    elseif strcmpi(Settings.paradigm, 'rest'),
        mats{iSession} = ROInets.run_individual_network_analysis(D,                          ...
                                                                 Settings,                   ...
                                                                 matsSaveFileName{iSession}, ...
                                                                 iSession);
    else
        error([mfilename ':BadParadigm'], ...
              'Unrecognised paradigm %s. \n', Settings.paradigm);
    end%if
end%for

% reformat results - correlationMats is a cell array of frequency bands
correlationMats = ROInets.reformat_results(mats, Settings);

% save current results: will write over later
% just in case of crash at group level
saveFileName = fullfile(outputDirectory, 'ROInetworks_correlation_mats.mat');
save(saveFileName, 'correlationMats');
clear mats

%% Subject-level analysis to average over sessions in a fixed-effects manner
correlationMats = ROInets.do_subject_level_glm(correlationMats, Settings);

%% Group-level analysis
% Find whole group means
if strcmpi(Settings.paradigm, 'rest'),
    correlationMats = ROInets.do_group_level_statistics(correlationMats, Settings);
end%if

% Perform group-level GLM
if ~isempty(Settings.GroupLevel),
    correlationMats = ROInets.do_group_level_glm(correlationMats, Settings);
end%if

%% save matrices
fprintf('\n%s: Saving Results. \n', mfilename);

% save collected results
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
