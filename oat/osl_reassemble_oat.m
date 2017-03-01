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
runcmd(['mkdir ' logdir]);

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
            runcmd(['cp ' subdir '/source_recon_mask.nii.gz ' masterDir]);
            
        end
        
        if (exist([subdir '/first_level_mask.nii.gz'],'file')) == 2
            % copy the source recon mask
            runcmd(['cp ' subdir '/first_level_mask.nii.gz ' masterDir]);
            
        end
    end
    
    % copy the files into the master oat directory
    runcmd(['cp ' subdir '/session' num2str(iSess) '_recon.mat ' masterDir]);
    
    runcmd(['cp ' subdir '/session' num2str(iSess) '_first_level.mat ' masterDir]);
    
    runcmd(['cp ' subdir '/concatMefsession' num2str(iSess) '_spm_meeg.mat ' masterDir]);
    
    runcmd(['cp ' subdir '/concatMefsession' num2str(iSess) '_spm_meeg.mat ' masterDir]);
    
    
    % copy the log file
    runcmd(['cp ' subdir '/log_sess_' num2str(iSess) '.txt ' logdir]);
    
    
    % delete the subdirectory
    runcmd(['rm -r ' subdir]);
    
    
    disp(['Adding fields to oat, session ' num2str(iSess) ' of ' num2str(nSess)]);
    %% add fields to the oat
    oat.source_recon.results_fnames{iSess} = ['session' num2str(iSess) '_recon'];
    oat.first_level.results_fnames{iSess}  = ['session' num2str(iSess) '_first_level'];
    
end % for iSess = 1:nSess

oat.to_do = [0 0 1 1];
oat = osl_save_oat(oat);

disp('Done.')