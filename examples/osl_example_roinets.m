%% Example script for ROInets-based analysis

osldir = getenv('OSLDIR');
data_dir = fullfile(osldir,'example_data','roinets_example');
output_directory = fullfile(osldir,'practical','roinets_demo');
mkdir(output_directory)

%%
% Copy the demo files and make a list of the D files
subjects = 1:10;
D_files = {};
mkdir(fullfile(osldir,'practical','roinets_demo'));
for j = 1:length(subjects)
    fname = sprintf('subject_%d',j);
    D = spm_eeg_load(fullfile(data_dir,fname));
    D_files{j} = D.copy(fullfile(output_directory,fname));
end

%% If files already exist
% subjects = 1:10;
% D_files = {};
% for j = 1:length(subjects)
%     D_files{j} = fullfile(data_dir,sprintf('subject_%d',j));
% end


%% 
% Load the parcellation and save a parcelflag variable
% Settings.spatialBasisSet = fullfile(output_directory,'parcelflag.nii.gz');
% p = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
% p.savenii(p.to_matrix(p.binarize),Settings.spatialBasisSet);
% Settings.gridStep = p.resolution;

%%
% Provide input parcellation
Settings = struct();
Settings.spatialBasisSet = fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
Settings.gridStep = p.resolution;
Settings = ROInets.check_inputs(Settings);
Settings.nSessions = numel(D_files);



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
