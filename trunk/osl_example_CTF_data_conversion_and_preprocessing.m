%% OSL CTF Data Conversion & Preprocessing Script Example 

% You will need the CTF 2-back data folder "OSL_example_CTF_raw_data"
% which includes: 
% and OSL version 1.1

% Henry Luckhoo 21.05.12

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

global OSLDIR;
    
%tilde='/home/mwoolrich';
tilde='/Users/woolrich';
osldir=[tilde '/homedir/matlab/osl2.0'];    

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

datadir=['/Users/woolrich/homedir/matlab/osl2_testdata_dir']; % directory where the data is
workingdir=[datadir '/ctf_fingertap_subject1_data_osl2']; % this is the directory where the SPM files will be stored in
cmd = ['mkdir ' workingdir]; unix(cmd); % make dir to put the results in
cd(workingdir);

clear ctf_files spm_files structural_files fidfiles

% list of CTF files
ctf_files={[workingdir '/JRH_MotorCon_20100429_01_FORMARK.ds']};

% set up a list of SPM MEEG object file names (we only have one here)
spm_files={[workingdir '/subject1.mat']};
   
ds_spm_files={[workingdir '/dsubject1.mat']};

% structural files in the same order as spm_files or .ds files:  
%structural_files = {[workingdir '/subject1_struct.nii']};      
structural_files = {[workingdir '/subject1_struct.nii']};      

% fiducial files in the same order as spm_files or .ds files:          
fidfiles={[workingdir '/subject1_Motorcon.pos']};
 
cleanup_files=0; % flag to indicate that you want to clean up files that are no longer needed

%%%%%%%%%%%%%%%%%%%%%%%%%
%% CTF data preprocessing
% CONVERT FROM .ds TO AN SPM MEEG OBJECT:


if(length(ctf_files)>0),
    % loads fif into SPM format (fif will typically have been maxfiltered and
    % downsamped (by maxfilter))
    % (and does detection of events)
    
    S2=[];

    for i=1:length(ctf_files), % iterates over subjects
        
        S2.fif_file=ctf_files{i};
        S2.spm_file=spm_files{i};       

        % The conversion to SPM will show a histogram of the event codes
        % and correspond to those listed below in the epoching section
        [D spm_files{i}] = osl_convert_script(S2);
    end;
end;

% Note that this spmfile is the output from batch_script1:
spm_files{1}

%%%%%%%%%%%%%%%%%%%
%% Downsample

for i=1:length(ctf_files), % iterates over subjects
        
    S2=[];
    S2.fsample_new=150;
    S2.D=spm_files{i};       

    % The conversion to SPM will show a histogram of the event codes
    % and correspond to those listed below in the epoching section
    D = spm_eeg_downsample(S2);
end;

%%%%%%%%%%%%%%%%%%%
%% DO REGISTRATION AND RUN FORWARD MODEL BASED ON STRUCTURAL SCANS
% Before running the beamformer we need to compute the forward model for
% each subject based on their structural scan.
% Make sure you check the results look reasonable!

for i=1:length(ds_spm_files),
    
        D=spm_eeg_load(ds_spm_files{i});
        
        try,D=rmfield(D,'inv');D.save;catch, end;
       
        
        %fidnew=D.fiducials;
        %fidnew.fid.label

        S2=[];
        
        S2.fid_label.nasion='nas';        
        S2.fid_label.lpa='lpa';
        S2.fid_label.rpa='rpa';
        
        S2.D = ds_spm_files{i};    % NOTE: requires .mat extension    

        
        S2.forward_meg='Single Shell';
        S2.mri=structural_files{i}; % set S2.sMRI=''; if there is no structural available        
        S2.useheadshape=1;
        S2.sMRI=S2.mri;
        
        D=osl_forward_model(S2);
end;

ca;
D=spm_eeg_load(ds_spm_files{1});
%spm_eeg_inv_checkmeshes(D);
spm_eeg_inv_checkdatareg(D);
%spm_eeg_inv_checkforward(D, 1);


