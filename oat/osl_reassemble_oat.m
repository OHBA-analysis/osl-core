function oat = osl_reassemble_oat(oat)
%
% Takes an oat computed in parallel, and reassembles it (for subject and
% group level)

osldir = '/home/gwallis/matlab/analysis/osl1.3.0';
addpath(osldir);
osl_startup(osldir);

%% load the 'master oat'
oat = osl_load_oat(oat);

% get the master oat directory
masterDir = oat.source_recon.dirname;
[rootDir,master] = fileparts(masterDir);

% make a folder to keep the logfiles from the source recon and first level,
% for debugging analyses
logdir = [masterDir '/logfiles'];
dos(['mkdir ' logdir]);

sprintf('Oat is in %s. \n\n Re-assembling the oat...',masterDir);

nSess = numel(oat.source_recon.sessions_to_do);
oat.source_recon.results_fnames = cell(1,nSess);
oat.first_level.results_fnames  = cell(1,nSess);
for iSess = 1:nSess
    
    disp(['Rearranging files,    session ' num2str(iSess) ' of ' num2str(nSess)]);
    %% File ops
    % construct the root of the subfolder for this session
    nam = [master '_sess' num2str(iSess) '.oat'];
    subdir = fullfile(rootDir,nam);
    
    if iSess == 1
        if (exist([subdir '/source_recon_mask.nii.gz'],'file')) == 2
            % copy the source recon mask
            cmmd = ['cp ' subdir '/source_recon_mask.nii.gz ' masterDir];
            [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
        end
        
        if (exist([subdir '/first_level_mask.nii.gz'],'file')) == 2
            % copy the source recon mask
            cmmd = ['cp ' subdir '/first_level_mask.nii.gz ' masterDir];
            [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
        end
    end
    
    % copy the files into the master oat directory
    cmmd = ['cp ' subdir '/session' num2str(iSess) '_recon.mat ' masterDir];
    [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
    cmmd = ['cp ' subdir '/session' num2str(iSess) '_first_level.mat ' masterDir];
    [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
    cmmd = ['cp ' subdir '/concatMefsession' num2str(iSess) '_spm_meeg.mat ' masterDir];
    [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
    cmmd = ['cp ' subdir '/concatMefsession' num2str(iSess) '_spm_meeg.mat ' masterDir];
    [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
    
    % copy the log file
    cmmd = ['cp ' subdir '/log_sess_' num2str(iSess) '.txt ' logdir];
    [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
    
    % delete the subdirectory
    cmmd = ['rm -r ' subdir];
    [s,~] = dos(cmmd); if s~=0; error(['Error with command: ' cmmd]); end
    
    disp(['Adding fields to oat, session ' num2str(iSess) ' of ' num2str(nSess)]);
    %% add fields to the oat
    oat.source_recon.results_fnames{iSess} = ['session' num2str(iSess) '_recon'];
    oat.first_level.results_fnames{iSess}  = ['session' num2str(iSess) '_first_level'];
    
end % for iSess = 1:nSess

oat.to_do = [0 0 1 1];
oat = osl_save_oat(oat);

disp('Done.')