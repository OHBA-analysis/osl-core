% TAKE IN - a directory
% And optionally specify the subject IDs? Session IDs? Or just the raw directory?
%
% Preprocessing consists of two stages
% 1. Import and coreg - where we have fif files, nifti files, and MEEG files
% 2. Standard preproc and beamforming - now we can deal with just a folder of MEEGs

% Form the input to the pipeline
subject_ids = {'DFPP04001','DFPP04002','DFPP04003','DFPP04004'};
conditions = {'RSC','RSO','SC1','SC2','MMN'};
S = [];
for j = 1:length(subject_ids)
	for k = 1:length(conditions)
		S(end+1).raw = sprintf('%s_%s.fif',subject_ids{j},conditions{k});
		S(end).nii = sprintf('%s_struc.nii.gz',subject_ids{j});
		S(end).meeg = sprintf('%s_%s.mat',subject_ids{j},conditions{k});

		% Now we create subfolders for each stage of the pipeline 
	end
end

% Remove nonexistant files, otherwise copy into working directory
failed = false(size(S));
for j = 1:length(S)
	if ~exist(S(j).raw)
		fprintf(2,'%s not found\n',S(j).raw);
		failed(j) = 1;
	end
end
S(failed) = [];
use_existing = true;


%% DOUBLE MAXFILTER AND IMPORT
import_dir = 'maxfilter';
if ~exist(import_dir,'dir')
    mkdir(import_dir)
end

D = cell(length(S),1);
for j = 1:length(S)

	if use_existing && exist(fullfile(import_dir,S(j).meeg),'file')
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

%% SINGLE MAXFILTER AND IMPORT
import_dir = 'no_maxfilter';
D = cell(length(S),1);
parfor j = 1:length(S)

	if use_existing && exist(fullfile(import_dir,S(j).meeg),'file')
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

%% COREGISTER
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
	report.save_figs(h,S(j).import,D{j}.fname);
	close(h);
end

%% NOW WE ARE DONE WITH IMPORTING AND WORKING WITH RAW FILES

%% Copy the inverse model from the no-sss ones to the sss

%% Copy and downsample
wd = 'no_maxfilter_analysis';
mkdir(wd);

s = study('no_maxfilter')
parfor j = 1:s.n
	if use_existing && exist(fullfile(wd,s.fnames{j}),'file')
		fprintf('Already downsampled %s\n',s.fnames{j});
		continue
	end

	D = s.read(j);
	fprintf('Downsampling %s\n',D.fname)
	spm_eeg_downsample(struct('D',D,'fsample_new',250,'prefix',[wd '/'])); % Note - downsampling cannot be done in-place using prefix='', it just fails
end

s = study(wd)
%% Do some filtering
D = osl_filter(D,[0.1 inf],'prefix',''); % Remove slow drift
D = osl_filter(D,-1*(50+[-2 2])); % Remove 50Hz with notch filter
D = osl_filter(D,-1*(2*100+[-2 2])); % Remove 100Hz with notch filter
D = osl_detect_artefacts(D);

S = struct;
S.artefact_channels = {'EOG1','ECG','EMG','EOG2'};
S.precompute_topos = false;
S.ident_func = @identify_artefactual_components_auto;
S.ident_params = struct('artefact_chans_corr_thresh',0.4); 
% You might also want to enable or disable other types of ICA component rejection. Finally, call |osl_africa|
D = osl_africa(D,S);
has_montage(D)



S = study('')




%% Detect artefacts
for j = 1:length(D)
	D{j} = osl_import(fullfile(S(j).dir,S(j).raw));
	h = report.events(D{j}); % Check the events have been read in
	report.save_figs(h,S(j).dir,D{j}.fname);
	close(h);
end


% And - what if we want to re-run each stage? Or could this just be designed as a script?
% Well - we would want to make it so that 

% BASIC USAGE
% - Takes in a pipeline file, that does D = pipeline(S), and S
% - where S is a struct array and each entry gets passed in turn to the pipeline
% - The pipeline expects to be working in an isolated subdirectory
% - Whether runs can be paused part way or whatnot is controlled by the contents of the pipeline file. 
%   Pipeline could decide to skip files if they already exist, or user could comment out bits and pieces of it
% - Pipeline returns a D if it exited correctly, or an MException if it did not (this is a way for errors to be passed up) or [] if the pipeline author was in a hurry
% SCENARIOS
% 
% 1. Rip through the data
%  - Construct list of files


DIRECTORY IN


%  - [output_dir,error_indices] = new_opt(S)
% 2. 




% 
% - INPUTS - struct array w list of fif files/DS folders, and structural files

% How do I provide inputs 



% The output is created IN THIS DIRECTORY
% The output consists of a folder for each FILE being processed that contains
% the intermediate STUFF
% And then a final output directory that stores the result of the pipeline
% That is, we pass in a function that takes in 


% The idea is - fif+structural goes in, and 



datadir = '/home/disk3/ajquinn/Projects/WakeHen/raw_data/';
procdir = '/home/disk3/ajquinn/Projects/WakeHen/preproc_data/';
subs = 1:19;

fif_files = cell(19*6,1);
conv_files = cell(19*6,1);
proc_files = cell(19*6,1);
epoched_files = cell(19*6,1);
beamformed_files = cell(19*6,1);
parc_files = cell(19*6,1);
closparc_files = cell(19*6,1);
struct_files = cell(19*6,1);

su = cell(19,1);se = cell(19,1);
for sub = subs
    for session = 1:6
        ind = (sub-1)*6 + session;
        su{ind} = sub;se{ind} = session';
        fif_files{ind} = [ datadir 'sub0' num2str(sprintf('%02d',sub)) '/run_0' num2str(session) '_sss.fif'];
        conv_files{ind} = [ procdir 'spm_sub' num2str(sprintf('%02d',sub)) '_run_0' num2str(session) '_sss.mat'];
        proc_files{ind} = [ procdir 'fdspm_sub' num2str(sprintf('%02d',sub)) '_run_0' num2str(session) '_sss.mat'];
        beamformed_files{ind} = [ procdir 'Bfdspm_sub' num2str(sprintf('%02d',sub)) '_run_0' num2str(session) '_sss.mat'];
        parc_files{ind} = [ procdir 'pBfdspm_sub' num2str(sprintf('%02d',sub)) '_run_0' num2str(session) '_sss.mat'];
        closparc_files{ind} = [ procdir 'cBfdspm_sub' num2str(sprintf('%02d',sub)) '_run_0' num2str(session) '_sss.mat'];
        struct_files{ind} = [ datadir 'sub0' num2str(sprintf('%02d',sub)) '/highres001.nii'];
    end
 end




