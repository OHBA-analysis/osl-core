% This practical we will work with a single subject's data from a blocked
% finger tapping experiment (CTF data courtesy of Matt Brookes), and
% perform a time-frequency analysis in the beta band using a time-wise GLM
% in sensor space.   
 
% Data can be downloaded from: 
% www.fmrib.ox.ac.uk/~woolrich/ctf_fingertap_subject1_data.tar.gz

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

OSLDIR = getenv('OSLDIR');
    
% set this to where you have downloaded OSL and the practical data:
practical_dir='/home/mwoolrich/Desktop'; 
osldir=[practical_dir '/osl1.3.1'];    

%practical_dir='/Users/woolrich';
%osldir=[practical_dir '/homedir/matlab/osl1.3.1'];    

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

% directory where the data is:
workingdir=[practical_dir '/ctf_fingertap_subject1_data']; % directory where the data is
%workingdir=[practical_dir '/homedir/matlab/osl_testdata_dir/ctf_fingertap_subject1_data']; % directory where the data is

cmd = ['mkdir ' workingdir]; unix(cmd); % make dir to put the results in

cd(workingdir);
% Set up the list of subjects and their structural scans for the analysi 
% Currently there is only 1 subject.
clear spm_files;

% set up a list of SPM MEEG object file names (we only have one here)
spm_files={[workingdir '/dsubject1.mat']};
 
cleanup_files=0; % flag to indicate that you want to clean up files that are no longer needed

%%%%%%%%%%%%%%%%%%%%
%% load the file in and take a look at the SPM object.

subnum = 1;             
D = spm_eeg_load(spm_files{subnum});

% look at the SPM object. Note that it is continuous data, with 216000 time
% points at 150Hz. 
D

%%%%%%%%%%%%%%%%%%%
%% Look at the data in oslview 
% Note that this is CTF data and so there is only one channel type (axial
% gradiometers).
% You will see that there are no obvious artefacts to mark.

oslview(D);

%%%%%%%%%%%%%%%%%%%
%% setup oat

oat=[];
oat.source_recon.D_continuous=spm_files;
oat.source_recon.conditions={'Undefined'};
oat.source_recon.freq_range=[13 30]; % frequency range in Hz
oat.source_recon.time_range=[300,32*30];
%oat.source_recon.time_range=[300,14*30];
oat.source_recon.method='none';
oat.source_recon.modalities={'MEGGRAD'};
oat.source_recon.work_in_pca_subspace=1;
oat.source_recon.pca_dim=270;

oat.source_recon.dirname=[workingdir '/subj1_results.oat'];

oat = osl_check_oat(oat);

%%%%%%%%%%%%%%%%%%%%%%%
%% run "source_recon" stage 
% Note that this is a sensor space analysis and so this just preps the data
% for the first level stage  

oat.to_do=[1 0 0 0];

oat = osl_run_oat(oat);

%%%%%%%%%%%%%%%
%% Establish regressor for continuous time GLM. 
% This should be setup to correspond to the time window
% over which the source recon is carried out.

D=spm_eeg_load(spm_files{1});

% This should be setup to correspond to the same time window
% as the full time window for D

x=zeros(length(D.time),5);

block_length=30; %s
block_order=[5 5 5 5 5 5 5 5 5 5 4 3 2 1 2 3 1 4 3 4 1 3 2 1 4 4 2 1 3 3 4 1 4 3 1 2 1 2 3 4 3 4 1 2 3 4 1 2];

% [Left, Right, Rest, Both, Rest_at_start]
% [  1     2      4     8     16 ]
% figure;plot(D.time,squeeze(D(1,:,:)))
% emacs ~/vols_data/From_Nottingham_with_Love/JRH_MotorCon_20100429_01_FORMARK.ds/MarkerFile.mrk 

tres=1/(D.fsample);
tim=1;
for tt=1:length(block_order),    
    x(tim:tim+block_length/tres-1,block_order(tt))=1;
    tim=tim+block_length/tres;
end;

figure;plot(D.time,x);

%%%%%%%%%
%% Run time-wise glm to do regression against known finger tapping regressors

oat=osl_load_oat(oat.source_recon.dirname);

oat.first_level.time_moving_av_win_size=1;
oat.first_level.design_matrix=x';
oat.first_level.contrast{1}=[-1 0 1 0 0]'; % rest-left
oat.first_level.contrast{2}=[0 -1 1 0 0]'; % rest-right
oat.first_level.contrast{3}=[0  0 1 -1 0]'; % rest-both
oat.first_level.contrast_name{1}='rest-left';
oat.first_level.contrast_name{2}='rest-right';
oat.first_level.contrast_name{3}='rest-both';
oat.first_level.bc=[0 0 0];
oat.first_level.diagnostic_cons_to_do=1:3;
oat.first_level.tf_freq_range=[13 30];
oat.first_level.tf_hilbert_freq_res=diff(oat.first_level.tf_freq_range);
oat.first_level.tf_method='hilbert';
oat.first_level.post_tf_downsample_factor=2; 

oat.first_level.doGLM=1;

oat.first_level.name=['subj1_first_level_ft'];

oat.to_do=[0 1 0 0];

oat = osl_run_oat(oat);

%% visualise using Fieldtrip
% note that this produces an interactive figure, with which you can:
% - draw around a set of sensors
% - click in the drawn box to produce a plot of the time series
% - on the time series plot you can draw a time window
% - and click in the window to create a topoplot averaged over that time
% window (which is itself interactive....!)

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGGRAD';
S2.first_level_contrast=1;
S2.cfg.colorbar='yes';
S2.cfg.interactive='yes';
%S2.cfg.zlim=[0 80];S2.view_cope=1;
S2.view_cope=0; % means we view the t-stats

% calculate t-stat using contrast of absolute value of parameter estimates
oat_stats_multiplotER(S2);

% try taking a look at the other contrasts and look at the lateralisation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at beta power time course at a sensor over motor cortex

stats=oat_load_results(oat,oat.first_level.results_fnames{1});
% chan_label='MRF54';chan_ind=find(strcmp(stats.D_sensor_data.chanlabels,chan_label));
[a chan_ind]=max(squeeze(stats.cope(:,1,S2.first_level_contrast)))

% re-run oat without doing GLM - this will output the time series used to
% fit the GLM to, i.e. the beta power (Hilbert envelope) time courses
oat.to_do=[0 1 0 0];
oat.first_level.doGLM=0;
oat = osl_run_oat(oat);

% do plot
stats2=oat_load_results(oat,oat.first_level.results_fnames{1});
figure;plot(stats2.glm_input_times,normalise(squeeze(stats2.glm_input_data(chan_ind,:))));
% compare to design matrix:
time_ind=intersect(find(D.time>=oat.source_recon.time_range(1)),find(D.time<=oat.source_recon.time_range(2)));
ho;plot(D.time(time_ind),x(time_ind,1),'r','LineWidth',2);
plot(D.time(time_ind),x(time_ind,2),'g','LineWidth',2);
plot(D.time(time_ind),x(time_ind,4),'k','LineWidth',2);
legend('data','left','right','both');
plot4paper('time(secs)','beta power');
