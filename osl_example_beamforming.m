% This practical will work with a single subject's data from an emotional
% faces experiment (Elekta Neuromag data). 
%
%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

global OSLDIR;

% set this to where you have downloaded OSL and the practical data:
practical_dir='/home/mwoolrich/Desktop';
osldir=[practical_dir '/osl2.0'];

practical_dir='/Users/woolrich';
osldir=[practical_dir '/Dropbox/osl2.0'];

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

% directory where the data is:
workingdir=[practical_dir '/homedir/vols_data/osl_testdata/osl2_testdata_dir/faces_subject1_data_osl2']; % directory where the data is

%workingdir = '/Users/andrew/Software/Matlab/osl_test_data/face_data_sub1_osl2_new';

cmd = ['mkdir ' workingdir]; if ~exist(workingdir, 'dir'), unix(cmd); end % make dir to put the results in

clear processed_files beamformed_files;

processed_files{1}=[workingdir '/dspm_meg1.mat'];

beamformed_files{1}=[workingdir '/bdspm_meg1.mat'];

%%%%%%%%
%% beamform

% Copy, bandpass filter & beamform
for session = 1:numel(processed_files)

    % Copy into beamforming directory:
    S           = [];
    S.D         = processed_files{session};
    S.outfile   = beamformed_files{session};
    spm_eeg_copy(S)
    
    % Band pass filter:
    S      = [];
    S.D    = beamformed_files{session};
    S.band = 'bandpass';
    S.freq = freq;
    osl_spmfun(@spm_eeg_filter,S)
    
    % Beamform:
    S                   = [];
    S.D                 = beamformed_files{session};
    S.modalities        = {'MEGGRAD'};
    S.mni_coords        = osl_mnimask2mnicoords(fullfile(OSLDIR,'/std_masks/MNI152_T1_8mm_brain.nii.gz'));
    S.timespan          = [0 Inf];
    S.pca_order         = 250;
    S.type              = 'Scalar';
    S.inverse_method    = 'beamform';
    S.prefix            = '';
    osl_inverse_model(S);

end
    
%%%%%%%%%%%%%%%%%%%%%%%
%% load in beamformer result

D=spm_eeg_load(beamformed_files{1});

% make sure that we switch from sensor space to source space:
D.montage('switch',2);
