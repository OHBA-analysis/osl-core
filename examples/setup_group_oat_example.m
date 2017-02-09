%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

% SET THE BELOW LINE TO THE OSL DIRECTORY:
tilde='/home/mwoolrich';
%tilde='/Users/woolrich';

osldir=[tilde '/homedir/matlab/osl1.5.0_beta'];
addpath(osldir);
osl_startup(osldir);


%%%%%%%%%%%%%%%%%%
%%

clear subjinitials;

subjinitials{1}='rmd'; subjdir{1}='/net/aton/data/perl_pilot/dutton_rebekah/090212'
subjinitials{2}='rnk'; subjdir{2}='/net/aton/data/perl_pilot/kanso_riam/090212'
subjinitials{3}='fnd'; subjdir{3}='/net/aton/data/perl_pilot/doran_fiona/090219'
subjinitials{4}='kab'; subjdir{4}='/net/aton/data/perl_pilot/burgess_kate/090226'

subjinitials{5}='zps'; subjdir{5}='/net/aton/data/perl_pilot/sherrell_zia/090319'
subjinitials{6}='fnf'; subjdir{6}='/net/aton/data/perl_pilot/fasanya_florence/090402' 
subjinitials{7}='bns'; subjdir{7}='/net/aton/data/perl_pilot/sunda_bobby/090416'
subjinitials{8}='mce'; subjdir{8}='/net/aton/data/perl_pilot/evans_matthew/090205'

subjinitials{9}='mdw'; subjdir{9}='/net/aton/data/perl_pilot/wilson_marcus/090226'
subjinitials{10}='fmc'; subjdir{10}='/net/aton/data/perl_pilot/cortella_fran/090312'
subjinitials{11}='dcd'; subjdir{11}='/net/aton/data/perl_pilot/dunlop_david/090312'
subjinitials{12}='ijs'; subjdir{12}='/net/aton/data/perl_pilot/simpson_ivor/090326'

subjinitials{13}='lpj'; subjdir{13}='/net/aton/data/perl_pilot/jones_lindsey/090326'
subjinitials{14}='fjw'; subjdir{14}='/net/aton/data/perl_pilot/woodward_felix/090402'

subjinitials{15}='eie'; subjdir{15}='/net/aton/data/perl_pilot/edstan_elisabeth/090129'
subjinitials{16}='lnj'; subjdir{16}='/net/aton/data/perl_pilot/jurkiewicz_louisa/090205'

subjinitials{17}='hjr'; subjdir{17}='/net/aton/data/perl_pilot/rose_hayley/090212'
subjinitials{18}='afd'; subjdir{18}='/net/aton/data/perl_pilot/dunham_april/090226'
subjinitials{19}='mic'; subjdir{19}='/net/aton/data/perl_pilot/chrisidu_marta/090226'
subjinitials{20}='lgd'; subjdir{20}='/net/aton/data/perl_pilot/dimitrov_lily/090423'

subjinitials{21}='hnt'; subjdir{21}='/net/aton/data/perl_pilot/townley_hannah/090423'
subjinitials{22}='tpf'; subjdir{22}='/net/aton/data/perl_pilot/furlong_tomas/090312'
subjinitials{23}='njm'; subjdir{23}='/net/aton/data/perl_pilot/moss_nicholas/090205'
subjinitials{24}='hif'; subjdir{24}='/net/aton/data/perl_pilot/fisher_harry/090305'

subjinitials{25}='mjs'; subjdir{25}='/net/aton/data/perl_pilot/sibley_marcus/090312'
subjinitials{26}='jnp'; subjdir{26}='/net/aton/data/perl_pilot/pinkney_justin/090319'
subjinitials{27}='ikl'; subjdir{27}='/net/aton/data/perl_pilot/vanderlowe_ilmo/090326'
subjinitials{28}='mdr'; subjdir{28}='/net/aton/data/perl_pilot/richard_mathew/090416'

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%% INITIALISE SETTINGS FOR THIS DATASET
% This specifies where our data is stored, what the data filenames are, and
% a few parameters for preprocessing such as how we will epoch the data.

basicdir=[tilde '/homedir/vols_data/nagels_project'];

fifdir=[basicdir '/fifs'];
workingdir=[basicdir '/spm_files'];
mkdir(workingdir);

%is=[1:length(subjinitials)];
is=15:28;

% Set up the list of subjects and their structural scans for the analysis. 
% Currently there is only 1 subject.
clear raw_fif_files structural_files spm_files_basenames;

% fif files:
raw_fif_files=[];
for i=1:length(is),    
    raw_fif_files{i}=[fifdir '/' subjinitials{is(i)} 'faceraw'];
end;

spm_files_basenames=[];
for i=1:length(is),    
    spm_files_basenames{i}=['spm8_meg' num2str(is(i))];    
end;

% structural files:
sfs={'F9'    'F11'    'F15'    'F16'    'F20'    'F22'    'F25'    'M3'    'M6'    'M7'    'M12'    'M14'    'M18'    'M19'    'F4'    'F7'    'F10'    'F13'    'F18'    'F24'    'F26'    'M2'    'M4'    'M8'    'M10'    'M13'   'M17'    'M21'    'F3'    'F6'    'F12'    'F14'    'F17'    'F19'    'F23'    'M5'    'M9'    'M11'    'M15'    'M20'    'M22'    'M23'};
sfs{21}=''; % structural is missing
sfs{14}=''; % structural is missing
sfs{6}=''; % structural is missing
for i=1:length(is),
    if ~strcmp(sfs{is(i)},'')
        structural_files{i}=[basicdir '/structurals/' sfs{is(i)} '/struct.nii'];
    else
        structural_files{i}='';
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%

testoutputdir=[tilde '/homedir/matlab/osl_testdata_dir/faces_group_data_osl1.5.1'];
mkdir(testoutputdir);

uber=[];
uber.do_sss_maxfilter=1;

uber.forward_meg='MEG Local Spheres';
uber.forward_meg='Single Shell';

uber.method='beamform';
uber.do_hmm=0;
uber.do_tf=1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%

opt=[];

tname=[osl2_version '_sss' num2str(uber.do_sss_maxfilter) '.opt'];

opt.dirname=[testoutputdir '/' tname];

opt.trigger_channel_mask='0000000000111111';

opt.raw_files=raw_fif_files;
opt.spm_files_basenames=spm_files_basenames;
opt.do_sss_maxfilter=uber.do_sss_maxfilter;
opt.do_remove_badchans_pre_sss=1;
opt.structural_files=structural_files;
opt.modalities={'MEGMAG';'MEGPLANAR'};

opt.downsample_freq=200;
opt.sss_movement_compensation=0;
opt.skip_all_maxfilter=0;

opt.epoch.time_range = [-1.2 1.2]; % epoch end in secs   

opt.epoch.trialdef(1).conditionlabel = 'Neutral face';
opt.epoch.trialdef(1).eventtype = 'STI101_down';
opt.epoch.trialdef(1).eventvalue = 1;
opt.epoch.trialdef(2).conditionlabel = 'Happy face';
opt.epoch.trialdef(2).eventtype = 'STI101_down';
opt.epoch.trialdef(2).eventvalue = 2;
opt.epoch.trialdef(3).conditionlabel = 'Fearful face';
opt.epoch.trialdef(3).eventtype = 'STI101_down';
opt.epoch.trialdef(3).eventvalue = 3;
opt.epoch.trialdef(4).conditionlabel = 'Motorbike';
opt.epoch.trialdef(4).eventtype = 'STI101_down';
opt.epoch.trialdef(4).eventvalue = 4;

opt.fid_label.nasion='Nasion';
opt.fid_label.lpa='LPA';
opt.fid_label.rpa='RPA';

opt.sessions_to_do=1:14;

if(0)
    % run OPT
    opt=osl_run_opt(opt);
else
    % load pre-run OPT
    optdir=['/home/mwoolrich/homedir/vols_data/nagels_project/spm_files/osl1.5.0norm_trafix_ica_sss1_st0.opt'];
    opt=osl_load_opt(optdir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
opt.sessions_to_do=1:14;
    
oat=[];

oat.source_recon.dirname=[testoutputdir '/beamform'];

oat.source_recon.D_continuous=[]; % do not have continuous files
oat.source_recon.D_epoched=opt.results.spm_files_epoched; % only epoched files
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};

        oat.source_recon.freq_range=[1 48]; % frequency range in Hz
        oat.source_recon.time_range=[-0.2 0.35];

oat.source_recon.method=uber.method;
oat.source_recon.gridstep=8; % in mm, using a lower resolution here than you would normally, for computational speed
oat.source_recon.mri=structural_files;
oat.source_recon.sessions_to_do=opt.sessions_to_do;
oat.source_recon.forward_meg=uber.forward_meg;
oat.source_recon.force_pca_dim=1;
oat.source_recon.pca_dim=55;
oat.source_recon.forward_meg='MEG Local Spheres';

if(uber.do_hmm)
    oat.source_recon.hmm_num_states=5;
    oat.source_recon.hmm_num_starts=4;
    oat.source_recon.hmm_pca_dim=40;
    oat.source_recon.hmm_av_class_occupancy=hmm_av_class_occupancy;
end;

Xsummary={};Xsummary{1}=[1 0 0 0];Xsummary{2}=[0 1 0 0];Xsummary{3}=[0 0 1 0];Xsummary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary=Xsummary;

% contrasts to be calculated:
oat.first_level.contrast={};
oat.first_level.contrast{1}=[3 0 0 0]'; % motorbikes
oat.first_level.contrast{2}=[0 1 1 1]'; % faces
oat.first_level.contrast{3}=[-3 1 1 1]'; % faces-motorbikes
oat.first_level.contrast_name{1}='motorbikes';
oat.first_level.contrast_name{2}='faces';
oat.first_level.contrast_name{3}='faces-motorbikes';
oat.first_level.cope_type='acope';
oat.first_level.tf_method='none'; % can be morlet or hilbert
oat.first_level.name=['first_level_' oat.first_level.tf_method];

% group-level stage

%oat = osl_load_oat(oat.source_recon.dirname,'first_level','sub_level','group_level'); 
oat.source_recon.sessions_to_do=opt.sessions_to_do;
oat.first_level.sessions_to_do=oat.source_recon.sessions_to_do;
oat.subject_level.subjects_to_do=oat.first_level.sessions_to_do;
oat.group_level.subjects_to_do=oat.first_level.sessions_to_do;

%oat.group_level.time_range=[0.145 0.155];
oat.group_level.space_average=0;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.mask_fname='';oat.group_level=rmfield(oat.group_level,'mask_fname');
oat.group_level.use_tstat=0;
oat.group_level.group_varcope_time_smooth_std=0;
oat.group_level.group_varcope_spatial_smooth_fwhm=40;
oat.group_level.name='group_level';

if(0)
    oat.group_level.group_design_matrix=zeros(2,28);
    oat.group_level.group_design_matrix(1,1:14)=1;
    oat.group_level.group_design_matrix(2,15:28)=1;
    oat.group_level.group_contrast=[];   
    oat.group_level.group_contrast{1}=[0 1]';  
    oat.group_level.group_contrast{2}=[1 0]';  
    oat.group_level.group_contrast{3}=[1 -1]';  
    oat.group_level.group_contrast{4}=[1 1]';
else
    oat.group_level.group_design_matrix=ones(1,14);
    oat.group_level.group_contrast=[];   
    oat.group_level.group_contrast{1}=[1]';     
end;

% run OAT
oat.to_do=[1 1 1 0];

oat = osl_check_oat(oat);



% oat = osl_run_oat(oat);



% load previously run oat:
oat = osl_load_oat(oat);









%%%%%%%%%%%%%%%%%%%
%% OAT using ROI from start with TF

% load in wholebrain erf oat from above and modify the settings
oat.source_recon.dirname=[testoutputdir '/beamform'];
oat.first_level.tf_method='none'; % can be morlet or hilbert
oat.first_level.name=['first_level_' oat.first_level.tf_method];
oat=osl_load_oat(oat);

% change settings
oat.to_do=[1 1 1 0];

oat.source_recon.dirname=[testoutputdir '/beamform_roi'];

oat.source_recon.time_range=[-0.2 0.4];
oat.source_recon.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];

oat.first_level.tf_method='hilbert'; % can be morlet or hilbert
oat.first_level.name=['first_level_' oat.first_level.tf_method];

oat.first_level.tf_num_freqs=14; % we are keeping this unusally low in the practical for the sake of speed
oat.first_level.tf_hilbert_freq_res=6;
oat.first_level.tf_freq_range=[2 40];
oat.first_level.time_range=[oat.source_recon.time_range(1)+0.05 oat.source_recon.time_range(2)-0.05];

oat.first_level.bc=[1 1 1];
oat.first_level.cope_type='cope';

oat.first_level.sessions_to_do=[1:14];


oat.group_level.time_range=oat.first_level.time_range;
oat.group_level.space_average=1;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
%oat.group_level.group_varcope_time_smooth_std=0.2;
%oat.group_level.group_varcope_spatial_smooth_fwhm=100;
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.use_tstat=0;
oat.group_level.name='group_level';

oat=osl_check_oat(oat);

%%

oat=osl_run_oat(oat);








%%%%%%%%%%%%%%%%%%%
%% OAT using sensor space
%% SETUP THE OAT:

oat=[];
oat.source_recon.D_continuous=[]; % do not have continuous files
oat.source_recon.D_epoched=opt.res.spm_files_epoched; % only epoched files
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.45];
oat.source_recon.method='none';
% load in wholebrain oat
tname=['faces_group_sensor_norm'];
oat.source_recon.dirname=[testoutputdir '/' tname];
oat.source_recon.sessions_to_do=[1:14];

% Xsummary is a parsimonious description of the design matrix.
% It contains values Xsummary{reg,cond}, where reg is a regressor no. and cond
% is a condition no. This will be used (by expanding the conditions over
% trials) to create the (num_regressors x num_trials) design matrix:
Xsummary={};
Xsummary{1}=[1 0 0 0];Xsummary{2}=[0 1 0 0];Xsummary{3}=[0 0 1 0];Xsummary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary=Xsummary;

% contrasts to be calculated:
oat.first_level.contrast={};
oat.first_level.contrast{1}=[3 0 0 0]'; % motorbikes
oat.first_level.contrast{2}=[0 1 1 1]'; % faces
oat.first_level.contrast{3}=[-3 1 1 1]'; % faces-motorbikes
%oat.first_level.contrast{4}=[0 1 0 0]'; % neutral
%oat.first_level.contrast{5}=[0 0 0 1]'; % fearful
%oat.first_level.contrast{6}=[0 -1 0 1]'; % fearful-neutral
oat.first_level.contrast_name={};
oat.first_level.contrast_name{1}='motorbikes';
oat.first_level.contrast_name{2}='faces';
oat.first_level.contrast_name{3}='faces-motorbikes';
%oat.first_level.contrast_name{4}='neutral';
%oat.first_level.contrast_name{5}='fearful';
%oat.first_level.contrast_name{6}='fearful-neutral';

oat.first_level.tf_method='hilbert'; % can be morlet or hilbert
oat.first_level.bc=[1 1 0];
oat.first_level.sessions_to_do=oat.source_recon.sessions_to_do;
oat.first_level.time_range=[-0.15 0.4];

%oat.first_level.tf_num_freqs=12; % we are keeping this unusally low in the practical for the sake of speed
%oat.first_level.tf_hilbert_freq_res=4;
% oat.first_level.tf_freq_range=[3 40];

oat.first_level.tf_num_freqs=1; % we are keeping this unusally low in the practical for the sake of speed
oat.first_level.tf_freq_range=[5 20];
oat.first_level.tf_hilbert_freq_res=diff(oat.first_level.tf_freq_range);

oat = osl_check_oat(oat);

%% run OAT

oat.do_plots=0;
oat.to_do=[1 1 1 0];
oat.group_level.group_varcope_time_smooth_std=0.2;

oat = osl_run_oat(oat);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE uploading_osl.txt in matlab dir



