
% Preprocessing consists of two stages
% 1. Import and coreg - where we have fif files, nifti files, and MEEG files
% 2. Standard preproc and beamforming - now we can deal with just a folder of MEEGs

%% Setup
% For the first stage, we need to specify a set of raw files (e.g. |.fif| or |.ds|).
% Our input consists of three pieces of information
%
% * The name of the file
% * The name of the corresponding structural NIFTI file
% * The name of the output MEEG object we want to write
% 
% You can write this in any way you like, but it's easy to specify this in a struct array.
% The code below creates such an array
subject_ids = {'DFPP04001','DFPP04002','DFPP04003','DFPP04004'};
conditions = {'RSC','RSO','SC1','SC2','MMN'};
S = [];
for j = 1:length(subject_ids)
    for k = 1:length(conditions)
        S(end+1).raw = sprintf('%s_%s.fif',subject_ids{j},conditions{k});
        S(end).nii = sprintf('structurals/%s_struc.nii.gz',subject_ids{j});
        S(end).meeg = sprintf('%s_%s.mat',subject_ids{j},conditions{k});
    end
end

%%
% In our case, some of the files are missing. So we need to go over the array
% and check whether the files are actually there or not.
failed = false(size(S));
for j = 1:length(S)
    if ~osl_util.isfile(S(j).raw) && ~osl_util.isfile(fullfile(import_dir,S(j).meeg))
        fprintf(2,'%s not found\n',S(j).raw);
        failed(j) = 1;
    end
end
S(failed) = [];

%%
% We are now ready to actually do the conversion.

%% Double maxfilter and import
% The first conversion approach we will take is to do the double maxfilter
%
% 1. Maxfilter 'nosss' to remove HPI signal
% 2. Detect bad channels
% 3. Maxfilter again using 'sss' with those channels removed
% 4. Import the 'sss' file 
%
% We will create a new directory for this, because we also want to try a
% pipeline without SSS. In general, it is a good idea to have only a single
% set of files in each folder. That is, it is convenient if a folder only
% contains one copy of the entire data set. So for example, we will have a
% folder for the maxfiltered data, and a folder for the non-maxfiltered data,
% rather than putting both data sets in the same folder. The reason for this
% should be clearer further in the pipeline when we start using the |study|
% object.
%
% First, create the output directory
import_dir = 'maxfilter_raw';
if ~isdir(import_dir)
    mkdir(import_dir)
end

%%
% The code below iterates over the files and runs the standard double maxfilter pipeline. There are
% two extra optional parts to it 
%
% * The |if use_existing| block checks if the final output file already exists, and skips the file entirely
% if it does. This means that you can re-run this cell and it will resume if the previous run was not completed
% for whatever reason (e.g. there was an error that stopped the processing). Similarly, if you want to re-run the import
% on a single file, you could delete that file, and re-run this cell, and it would only re-import the file that you deleted.
% * We use the 'report.' functions to generate records and summary information. Unlike in |opt|, this is performed manually 
% and explicitly. This also means you can regenerate the report figures at any time. The |report| functions return an array of 
% figure handles. You can pass them to |report.save_figs| to save the figures to disk.
%
D = cell(length(S),1);
use_existing = true;
for j = 1:length(S)

    if use_existing && osl_util.isfile(fullfile(import_dir,S(j).meeg))
        fprintf('Already imported %s\n',S(j).meeg);
        D{j} = spm_eeg_load(fullfile(import_dir,S(j).meeg));
        continue
    end

    nosss_fif = osl_maxfilter(S(j).raw,'nosss','verbose',true,'fif_out',fullfile(import_dir,['nosss_',S(j).raw]));
    D_temp = osl_import(nosss_fif);
    D_temp = osl_detect_artefacts(D_temp,'badtimes',false);
    h = report.bad_channels(D_temp);
    report.save_figs(h,import_dir,D_temp.fname);
    close(h);
    
    [sss_fif,bad_segments,headpos_file] = osl_maxfilter(S(j).raw,'sss','badchannels',D_temp.chanlabels(D_temp.badchannels),'verbose',true,'fif_out',fullfile(import_dir,['sss_',S(j).raw]));
    h = report.headpos(headpos_file); % Check the events have been read in
    report.save_figs(h,import_dir,S(j).meeg);
    
    [ds_fif,bad_segments] = osl_maxfilter(sss_fif,'nosss','verbose',true,'fif_out',fullfile(import_dir,['ds_',S(j).raw]));
    D{j} = osl_import(ds_fif,'bad_segments',bad_segments,'outfile',fullfile(import_dir,S(j).meeg));
    h = report.events(D{j}); % Check the events have been read in
    report.save_figs(h,import_dir,D{j}.fname);
    close(h);
    delete(D_temp); % Delete the temporary nosss MEEG
    delete(nosss_fif); % Delete the nosss FIF file
end

%% Perform coregistration
% Coregistration is similar to before - you specify a struct with the coregistration settings and then
% run |osl_headmodel|. This will add a |.inv| field to the MEEG object. 
use_existing = true;
for j = 1:length(S)

    if use_existing & isfield(D{j},'inv')
        fprintf('Already coregistered %s\n',S(j).meeg);
        continue
    end

    coreg_settings = struct;
    coreg_settings.D = D{j}.fullfile;
    coreg_settings.mri = S(j).nii;
    coreg_settings.useheadshape = true;
    coreg_settings.forward_meg = 'MEG Local Spheres';
    coreg_settings.use_rhino = true;
    coreg_settings.fid.label.nasion='Nasion';
    coreg_settings.fid.label.lpa='LPA';
    coreg_settings.fid.label.rpa='RPA';
    D{j} = osl_headmodel(coreg_settings);
    h = report.coreg(D{j});
    report.save_figs(h,import_dir,D{j}.fname);
    close(h);
end

%% SINGLE MAXFILTER AND IMPORT
% If we wanted to test not running SSS, we would instead run this cell, which only uses |nosss|, followed by the coregistration cell above. Of course, in normal usage, because
% the first stage of double maxfilter is to run |nosss|, we could just use the |nosss*.fif| files that are generated anyway as part
% of the double maxfilter pipeline. This cell would be used if you didn't want to run SSS at all.
import_dir = 'no_maxfilter';
use_existing = true;
D = cell(length(S),1);
parfor j = 1:length(S)

    if use_existing && osl_util.isfile(fullfile(import_dir,S(j).meeg))
        fprintf('Already imported %s\n',S(j).meeg);
        D{j} = spm_eeg_load(fullfile(import_dir,S(j).meeg));
        continue
    end

    nosss_fif = osl_maxfilter(S(j).raw,'nosss','verbose',true,'fif_out',fullfile(import_dir,['nosss_',S(j).raw]));
    D{j} = osl_import(nosss_fif,'outfile',fullfile(import_dir,S(j).meeg));
    h = report.events(D{j}); % Check the events have been read in
    report.save_figs(h,import_dir,D{j}.fname);
    close(h);
end

%% COPY COREGISTRATION
% We now have two folders with a complete set of MEEG objects. However, only one of them has been coregistered. Coregistration is
% stochastic and running it again would result in a slightly different coregistration, which would be a confound if we are interested in
% comparing two different pipelines applied to the same data. So the solution is to copy the coregistration from the first set to the second set
coregistered_dir = 'maxfilter_raw';
new_dataset_dir = 'no_maxfilter';
for j = 1:length(S)
    D_coreg = spm_eeg_load(fullfile(coregistered_dir,S(j).meeg)); % Load the coregistered MEEG
    D_todo = spm_eeg_load(fullfile(new_dataset_dir,S(j).meeg)); % Load the target MEEG
    D_todo.inv = D_coreg.inv; % Copy the coregistration/head model
    D_todo.save; % Commmit to disk
end

%% Clear workspace
% We are done with the struct array of filenames. Now we will use the |study| object to iterate over
% folders containing a single set of MEEG objects - this has the advantage that we no longer need to
% explicitly specify the filenames. To demonstrate this, we can clear the workspace to show that there
% are no remaining dependencies from the previous importing/coreg stage.
clear all

%% Copy and downsample
% The SSS data should be downsampled to make processing it more tractable
use_existing = true;
wd = 'maxfilter';
mkdir(wd);
s = study('maxfilter_raw')
parfor j = 1:s.n
    if use_existing && osl_util.isfile(fullfile(wd,s.fnames{j}))
        fprintf('Already downsampled %s\n',s.fnames{j});
        continue
    end

    D = s.read(j);
    fprintf('Downsampling %s\n',D.fname)
    spm_eeg_downsample(struct('D',D,'fsample_new',250,'prefix',[wd '/'])); % Note - downsampling cannot be done in-place using prefix='', it just fails
end

%% Do processing
wd = 'maxfilter_processed';
mkdir(wd);
s = study('maxfilter')
for j = 1:s.n
    D = s.read(j);

    %% Do some filtering
    D = osl_filter(D,[0.1 inf],'prefix',fullfile(pwd,[wd '/'])); % Remove slow drift. Note the use of fullfile(pwd,...) to specify an absolute path, as required by |osl_inverse_model|
    D = osl_filter(D,-1*(50+[-2 2]),'prefix',''); % Remove 50Hz with notch filter
    D = osl_filter(D,[1 45],'prefix','');% Seems to benefit
    D = osl_detect_artefacts(D);
    D = osl_africa(D,'precompute_topos',false,'used_maxfilter',true);

    %% Do epoching
    trialdef = struct('conditionlabel',{},'eventtype',{},'eventvalue',{});
    if strfind(D.fname,'MMN')
        % Want to label the eventtype-eventvalue pair as a particular condition
        trialdef(end+1) = struct('conditionlabel' , 'Dur25_500'        , 'eventtype' , 'STI101_down' , 'eventvalue' , 1);
        trialdef(end+1) = struct('conditionlabel' , 'Freq_450'         , 'eventtype' , 'STI101_down' , 'eventvalue' , 2);
        trialdef(end+1) = struct('conditionlabel' , 'Freq_550'         , 'eventtype' , 'STI101_down' , 'eventvalue' , 3);
        trialdef(end+1) = struct('conditionlabel' , 'Gap_500'          , 'eventtype' , 'STI101_down' , 'eventvalue' , 4);
        trialdef(end+1) = struct('conditionlabel' , 'intensity_10_500' , 'eventtype' , 'STI101_down' , 'eventvalue' , 5);
        trialdef(end+1) = struct('conditionlabel' , 'intensity10_500'  , 'eventtype' , 'STI101_down' , 'eventvalue' , 6);
        trialdef(end+1) = struct('conditionlabel' , 'LocLeft_500'      , 'eventtype' , 'STI101_down' , 'eventvalue' , 7);
        trialdef(end+1) = struct('conditionlabel' , 'LocRight_500'     , 'eventtype' , 'STI101_down' , 'eventvalue' , 8);
        trialdef(end+1) = struct('conditionlabel' , 'typical'          , 'eventtype' , 'STI101_down' , 'eventvalue' , 11);
    elseif strfind(D.fname, '_SC')
        trialdef(end+1) = struct('conditionlabel' , 'Novel'            , 'eventtype' , 'STI101_down' , 'eventvalue' , 101);
        trialdef(end+1) = struct('conditionlabel' , 'Repeat'           , 'eventtype' , 'STI101_down' , 'eventvalue' , 102);
        trialdef(end+1) = struct('conditionlabel' , 'Moon'             , 'eventtype' , 'STI101_down' , 'eventvalue' , 104);
    end

    if isempty(strfind(D.fname,'_RS'))
        D = osl_epoch(struct('prefix','','D',D,'trialdef',trialdef,'timewin',[-200 1000]))
    end

    p = parcellation('dk_cortical.nii.gz');
    D = osl_inverse_model(D,p.template_coordinates);
end










has_montage(D)



%% Detect artefacts
for j = 1:length(D)
    D{j} = osl_import(fullfile(S(j).dir,S(j).raw));
    h = report.events(D{j}); % Check the events have been read in
    report.save_figs(h,S(j).dir,D{j}.fname);
    close(h);
end




