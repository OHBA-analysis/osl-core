%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

OSLDIR = getenv('OSLDIR');
    
% set this to where you have downloaded OSL and the practical data:
practical_dir='/home/mwoolrich/Desktop'; 
osldir=[practical_dir '/osl1.3.1'];    

practical_dir='/Users/woolrich';
osldir=[practical_dir '/homedir/matlab/osl1.4.0'];    

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

% directory where the data is:
workingdir=[practical_dir '/homedir/matlab/osl_testdata_dir/faces_group_data_osl1.4.0']; % directory where the data is

oatdir=[workingdir '/none.oat'];

%%%%%%%%%%%%%%%%%%
%% load pre-run OAT analysis and run group-level stage for SENSOR space
%% analysis

% Load OAT analysis for which the first 2 stages (source recon and
% first-level GLM) have already been run.
% Note that the 1st level contrasts run are:
% S2.contrast{1}=[3 0 0 0]'; % motorbikes
% S2.contrast{2}=[0 1 1 1]'; % faces
% S2.contrast{3}=[-3 1 1 1]'; % faces-motorbikes

oat = osl_load_oat(oatdir,'first_level_none','sub_level','group_level'); 

%% run OAT

% oat.first_level=rmfield(oat.first_level,'diagnostic_cons_to_do');oat.first_level=rmfield(oat.first_level,'diagnostic_modality_to_do')
%  oat.group_level.report=rmfield(oat.group_level.report,'group_cons_to_do'); oat.group_level=rmfield(oat.group_level,'diagnostic_group_cons_to_do'); oat.group_level=rmfield(oat.group_level,'diagnostic_first_level_cons_to_do'); oat.group_level=rmfield(oat.group_level,'diagnostic_modality_to_do')


oat.group_level.group_design_matrix=[];
oat.group_level.group_design_matrix(1,:)=[1 1 1 1 1 1 1 1 1 1 1 1 1 1];
oat.group_level.group_contrast=[];
oat.group_level.group_contrast{1}=[1]'; % group average

if(0)
    Num_subjects=size(oat.group_level.group_design_matrix,2);
    
    oat.group_level.group_design_matrix=ones(2,Num_subjects); % first regressor models group average
    %oat.group_level.group_design_matrix(2,:)=demean([randn(Num_subjects,1)]); % behavioural measures for each subject 
    oat.group_level.group_design_matrix(1,1:6)=0;
    oat.group_level.group_design_matrix(2,7:end)=0;
    
    oat.group_level.group_contrast=[];
    oat.group_level.group_contrast{1}=[1 0]'; % group average
    oat.group_level.group_contrast{2}=[0 1]'; % behavioural measure effect
    oat.group_level.group_contrast{3}=[1 1]'; % average effect!
end;
%oat.group_level=rmfield(oat.group_level,'diagnostic_cons_to_do');

oat.group_level.group_varcope_time_smooth_std=0.1;
oat.group_level.space_average=0;
oat.do_plots=0;
oat.to_do=[0 0 0 1];

oat = osl_check_oat(oat);

oat = osl_run_oat(oat);

%%

oat.group_level.report.first_level_cons_to_do=[3 1 2];
oat.group_level.report.modality_to_do='MEGMAG';
report = oat_group_level_stats_report(oat,oat.group_level.results_fnames);

%% visualise group result using Fieldtrip
% note that this produces an interactive figure, with which you can:
% - draw around a set of sensors
% - click in the drawn box to produce a plot of the time series
% - on the time series plot you can draw a time window
% - and click in the window to create a topoplot averaged over that time
% window (which is itself interactive....!)

S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.modality='MEGPLANAR';
S2.first_level_contrast=3;
S2.group_level_contrast=1;
S2.cfg.colorbar='yes';
S2.view_cope=0;

% calculate t-stat using contrast of absolute value of parameter estimates
[cfg, data]=oat_stats_multiplotER(S2);

%% do cluster stats in sensor space

S=[];
S.oat=oat;
S.cluster_stats_thresh=4;
S.cluster_stats_nperms=500;
S.first_level_copes_to_do=[3];
S.group_level_copes_to_do=[1];
S.group_varcope_time_smooth_std=oat.group_level.group_varcope_time_smooth_std;
S.modality='MEGMAG';
[corrp ts] = oat_cluster_perm_sensor_tf(S);

gcon=1;
con=3;
figure;plot(corrp(:,:,con,:,gcon)')
figure;plot(ts(:,:,con,:,gcon)')
figure;plot(ts(:,:,con,:,gcon)'>S.cluster_stats_thresh)

%% plot corrected P-values

S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.modality=S.modality; % MUST be set to MEGMAG to look at oat_cluster_perm_sensor_tf as even as MEGPLANAR are viewed this way
S2.data=corrp(:,:,con,:,gcon);
S2.first_level_contrast=con;
%S2.data=ts(:,:,con,:,gcon)>S.cluster_stats_thresh;
S2.cfg.colorbar='yes';
[cfg, data]=oat_stats_multiplotER(S2);











%%%%%%%%%%%%%%%%%%%
%% Whole of sensor-space Time-Freq analysis

oat = osl_load_oat(oatdir,'first_level_hilbert','sub_level','group_level'); 

oat.to_do=[0 0 0 1];

oat.group_level.time_range=oat.first_level.time_range;
oat.group_level.space_average=1;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.use_tstat=0;
oat.group_level.name='group_level_sensorav';
oat.group_level.store_lower_level_copes=1;

oat=osl_check_oat(oat);

oat=osl_run_oat(oat);

%% Plot the group results

% load GLM result
stats=oat_load_results(oat,oat.group_level.results_fnames);

figure;
freqbin=1;
cons=1:3;cols={'r','g','b'};
for c=1:length(cons),
    con=cons(c); plot(stats.times, squeeze(mean(stats.cope(1,:,con,freqbin)./stats.stdcope(1,:,con,freqbin),4)),cols{c},'LineWidth',2); hold on;
end;
legend('Motorbikes','Faces','Faces-Motorbikes');
xlabel('time (s)'); ylabel('1-tailed t-stat'); title([num2str(stats.frequencies(freqbin)) ' Hz']);

%% Do 2D cluster permutation stats

% load GLM result
stats=oat_load_results(oat,oat.group_level.results_fnames);

figure;cons=1:3;
cols={'r','g','b'};
freqbin=1;
for c=1:length(cons),
    contrast=cons(c);
    cluster_forming_threshold=2.8;
    num_perms=500;
    tmp=stats.lower_level_copes{contrast};
    varcope_time_smooth_std=oat.group_level.group_varcope_time_smooth_std;
    tres=stats.times(2)-stats.times(1);
    
    corrp = osl_clustertf(permute(tmp,[2 4 3 1]),cluster_forming_threshold,num_perms,26,varcope_time_smooth_std,tres);
    
    con=cons(c); 
    subplot(121);plot(stats.times, squeeze(corrp),cols{c},'LineWidth',2); hold on; a=axis; a(3)=0.8;a(4)=1.2;axis(a);
    xlabel('time (s)'); ylabel('corrected P-value'); title([num2str(stats.frequencies(freqbin)) ' Hz']);
    subplot(122);plot(stats.times, squeeze(mean(stats.cope(1,:,con,freqbin)./stats.stdcope(1,:,con,freqbin),4)),cols{c},'LineWidth',2); hold on;
    xlabel('time (s)'); ylabel('1-tailed t-stat'); title([num2str(stats.frequencies(freqbin)) ' Hz']);

end;
legend('Motorbikes','Faces','Faces-Motorbikes');

