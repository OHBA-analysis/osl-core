%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

% SET THE BELOW LINE TO THE OSL DIRECTORY:
tilde='/home/mwoolrich';
%tilde='/Users/woolrich';

osldir=[tilde '/homedir/matlab/osl1.5.0_beta'];
addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

is=15:28;

basicdir=[tilde '/homedir/vols_data/nagels_project'];

%opt=osl_load_opt('/home/disk3/mwoolrich/vols_data/nagels_project/spm_files/osl1.4.0_sss1_st0.opt','opt-27-Jan-2014');

%opt=osl_load_opt('/home/mwoolrich/homedir/vols_data/nagels_project/spm_files/osl1.4.0norm_trafix_ica_sss1_st0.opt','opt-23-Jul-2014');
opt=osl_load_opt('/home/mwoolrich/homedir/vols_data/nagels_project/spm_files/osl1.5.0norm_trafix_ica_sss1_st0.opt');

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
    
    spm_epoched_files{i}=opt.results.spm_files_epoched{is(i)};
end;

workingdir=['/home/mwoolrich/homedir/matlab/osl_testdata_dir/faces_group_data_' osl2_version];
mkdir(workingdir);

%%

oat=[];

oat.source_recon.time_range=[-0.15 0.3];

oat.source_recon.work_in_pca_subspace=0;
oat.source_recon.D_continuous=[]; % do not have continuous files
oat.source_recon.D_epoched=spm_epoched_files; % only epoched files
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[1 48]; % frequency range in Hz

oat.source_recon.gridstep=9; % in mm, using a lower resolution here than you would normally, for computational speed
oat.source_recon.mri=structural_files;
oat.source_recon.forward_meg='Single Shell';
oat.source_recon.neuromag_planar_baseline_correction='none';

Xsummary={};Xsummary{1}=[1 0 0 0];Xsummary{2}=[0 1 0 0];Xsummary{3}=[0 0 1 0];Xsummary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary=Xsummary;

% contrasts to be calculated:
oat.first_level.contrast={};
oat.first_level.contrast{1}=[3 0 0 0]'; % motorbikes
oat.first_level.contrast{2}=[0 1 1 1]'; % faces
oat.first_level.contrast{3}=[-3 1 1 1]'; % faces-motorbikes
oat.first_level.contrast_name={};
oat.first_level.contrast_name{1}='motorbikes';
oat.first_level.contrast_name{2}='faces';
oat.first_level.contrast_name{3}='faces-motorbikes';
oat.first_level.cope_type='acope';

oat = osl_check_oat(oat);

% group-level stage
oat.group_level.space_average=0;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.mask_fname='';oat.group_level=rmfield(oat.group_level,'mask_fname');
oat.group_level.use_tstat=0;
oat.group_level.group_varcope_time_smooth_std=0;
oat.group_level.group_varcope_spatial_smooth_fwhm=1000;
oat.group_level.first_level_contrasts_to_do=1:3;
    
oat.group_level.group_design_matrix=ones(1,length(oat.source_recon.D_epoched));
oat.group_level.group_contrast=[];   
oat.group_level.group_contrast{1}=[1]'; 
    
oat.group_level.group_contrast_name=[];
oat.group_level.group_contrast_name{1}='groupmean';
    
oat.group_level.report.first_level_cons_to_do=3;
oat.group_level.store_lower_level_copes=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensor space, whole brain, ERF

oat.source_recon.method='none';
oat.first_level.tf_method='none';

oat.first_level.name=['first_level_' oat.first_level.tf_method];
oat.source_recon.dirname=[workingdir '/' oat.source_recon.method];

oat.to_do=[1 1 0 0];

oat = osl_check_oat(oat);

% oat=oat_consolidate_results(oat);
oat = osl_run_oat(oat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensor space, whole brain, TF

oat.source_recon.method='none';

oat.first_level.tf_method='hilbert';
oat.first_level.tf_hilbert_freq_res=8;
oat.first_level.time_range=[-0.1 0.25];
oat.first_level.tf_num_freqs=10;
oat.first_level.tf_freq_range=[4 44];

oat.first_level.name=['first_level_' oat.first_level.tf_method];
oat.source_recon.dirname=[workingdir '/' oat.source_recon.method];

oat.to_do=[0 1 1 0];

oat = osl_check_oat(oat);

% oat=oat_consolidate_results(oat);
oat = osl_run_oat(oat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamform, whole brain, ERF

oat.source_recon.method='beamform';
oat.first_level.tf_method='none';

oat.first_level.name=['first_level_' oat.first_level.tf_method];
oat.source_recon.dirname=[workingdir '/' oat.source_recon.method];

oat.to_do=[1 1 1 0];

oat = osl_check_oat(oat);

% oat=oat_consolidate_results(oat);
oat = osl_run_oat(oat);

%%
if(0)
    S2=[];
    S2.oat=oat;
    S2.stats_fname=oat.first_level.results_fnames{3};
    S2.first_level_contrasts=[3]; % list of first level contrasts to output
    S2.resamp_gridstep=oat.source_recon.gridstep;

    [statsdir,times]=oat_save_nii_stats(S2);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamform, mask (from first level), ERF
if(0)
    oat.source_recon.method='beamform';
    oat.first_level.tf_method='none';

    oat.first_level.name=['first_level_' oat.first_level.tf_method '_mask'];
    oat.source_recon.dirname=[workingdir '/' oat.source_recon.method];

    oat.first_level.mask_fname=[osldir '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];

    oat.to_do=[0 1 1 0];

    oat = osl_check_oat(oat);

    % oat=oat_consolidate_results(oat);
    oat = osl_run_oat(oat);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamform, mask (from first level), TF

oat.source_recon.method='beamform';

oat.first_level.tf_method='hilbert';
oat.first_level.tf_hilbert_freq_res=8;
oat.first_level.time_range=[-0.1 0.25];
oat.first_level.tf_num_freqs=10;
oat.first_level.tf_freq_range=[4 44];

oat.first_level.name=['first_level_' oat.first_level.tf_method '_mask'];
oat.source_recon.dirname=[workingdir '/' oat.source_recon.method];

oat.first_level.mask_fname=[osldir '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];

oat.to_do=[0 1 1 0];

oat = osl_check_oat(oat);

% oat=oat_consolidate_results(oat);
oat = osl_run_oat(oat);
















