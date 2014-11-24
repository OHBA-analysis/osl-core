function [results_fnames, source_recon_results ]=osl_run_source_recon_inverse(oat)

% results_fnames=osl_run_source_recon_beamform(oat)
%
% takes in an OAT, which needs to be setup by calling oat=osl_setup_oat(S), struct 
% and runs beamformer
% 
% This function should normally be called using osl_run_oat(oat);
%
% MW 2011

global OSLDIR;

source_recon=oat.source_recon;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup mni coords using std space mask
if(isfield(source_recon,'mask_fname'))

    if ~isempty(source_recon.gridstep)
        mask_fname=osl_resample_nii(source_recon.mask_fname, [source_recon.mask_fname '_' num2str(source_recon.gridstep) 'mm'], source_recon.gridstep, 'nearestneighbour', [OSLDIR '/std_masks/MNI152_T1_' num2str(source_recon.gridstep) 'mm_brain']);
    else
        mask_fname = source_recon.mask_fname;
    end

    disp(['Using mask ' mask_fname]);

    % setup mni coords using std space mask
    [ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname);
    
    str=['Using dipoles at MNI coordinates ' ];

    for vox=1:size(mni_coords,1),
        str=[str ', [' num2str(mni_coords(vox,:)) ']'];
    end;

    disp(str);
elseif(isfield(source_recon,'single_mni_coord'))
    disp(['Using MNI coordinate ' num2str(source_recon.single_mni_coord)]);
    mni_coords=source_recon.single_mni_coord;
    
    mask_fname='';
elseif(isfield(source_recon,'mni_coords'))    
    mni_coords=source_recon.mni_coords;
    
    str=['Using dipoles at MNI coordinates ' ];

    for vox=1:size(mni_coords,1),
        str=[str ', [' num2str(mni_coords(vox,:)) ']'];
    end;

    disp(str);
    
    mask_fname='';    
else
    mask_fname=[OSLDIR '/std_masks/MNI152_T1_' num2str(source_recon.gridstep) 'mm_brain'];
    disp(['Using whole brain']);
    
    % setup mni coords usource_recong std space mask
    [ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname);    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save mask in oat dir
if(~strcmp(mask_fname,''))
   mask=read_avw(mask_fname); 
   save_avw(mask,[oat.source_recon.dirname '/source_recon_mask'], 'f', [source_recon.gridstep,source_recon.gridstep,source_recon.gridstep,1]); 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set first level diagnostic report up    
report_dir=[oat.results.plotsdir '/' oat.results.date '_source_recon'];
source_recon_report=osl_report_setup(report_dir,['Source recon (epoched)']);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp(['Using ' source_recon.method ' for source recontruction']);  % changed dy DM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over sessions
for sessi_todo=1:length(source_recon.sessions_to_do),   
        
    sessi = source_recon.sessions_to_do(sessi_todo);
 
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OAT SOURCE RECON ON SESS = ' num2str(sessi) '  %%%%%%%%%%%%%%%%%%%%%%%'])
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set session specific diagnostic report up    
    report_dir = [source_recon_report.dir '/sess' num2str(sessi)];
    report     = osl_report_setup(report_dir,['Session ' num2str(sessi)]);       

    source_recon.session_name = oat.source_recon.session_names{sessi}; % set in osl_check_oat. Changed by GC to allow different naming conventions Jan 2014 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% setup source_recon_sess for this session
    source_recon_sess          = source_recon;
    source_recon_sess.do_plots = oat.do_plots;

    if ~isempty(source_recon.D_continuous),        
        [p, fname, e] = fileparts(source_recon.D_continuous{sessi});       
        source_recon_sess.D_continuous = [p '/' fname '.mat'];       
        disp('Using continuous data as input');        
    else
        source_recon_sess.D_continuous = [];
    end;
    
    if ~isempty(source_recon.D_epoched),
        [p, fname, e] = fileparts(source_recon.D_epoched{sessi});       
        source_recon_sess.D_epoched = [p '/' fname '.mat'];       
        disp('Using epoched data as input');
    else
        source_recon_sess.D_epoched = [];
    end;
    
    source_recon_sess.mri       = source_recon.mri{sessi};
    
    if length(source_recon.pca_dim)>1,
        source_recon_sess.pca_dim=source_recon.pca_dim(sessi);
    else
        source_recon_sess.pca_dim=source_recon.pca_dim;        
    end;
    
    source_recon_sess.mni_coords = mni_coords;
    
    if(isfield(source_recon,'hmm_block'))
        source_recon_sess.hmm_block = source_recon.hmm_block{sessi};
    end;
    
    source_recon_sess.report = report;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Prepare data (bandpass filter, epoch, normalise sensors, HMM, map into PCA)
    source_recon_results=[];
    [source_recon_sess source_recon_results report] = oat_prepare_source_recon(source_recon_sess, source_recon_results, report);
    D=source_recon_results.BF.data.D;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DO CO-REGISTRATION AND HEAD MODEL
    % Need to re-do it now following any montaging done in osl_prepare_oat_source_recon

    disp('Co-registration and setting up forward model...');

    S2 = [];
    S2.D            = [D.path '/' D.fname]; % requires .mat extension
    S2.mri          = source_recon.mri{sessi};
    S2.useheadshape = source_recon.useheadshape;
    S2.forward_meg  = source_recon.forward_meg;
    S2.fid_label    = source_recon.fid_label;

    D = osl_forward_model(S2);

    if(oat.do_plots)
        %spm_eeg_inv_checkmeshes(D);
        spm_eeg_inv_checkdatareg(D);
        %spm_eeg_inv_checkforward(D, 1);
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run inverse
    
    S2=[];
    S2.D=D;
    S2.mni_coords=mni_coords;
    S2.timespan=source_recon.time_range;
    S2.type=source_recon.type;
    S2.pca_order=source_recon.pca_dim;
    S2.modalities=source_recon.modalities;
    S2.fuse='all';
    S2.conditions=source_recon.conditions;
    S2.inverse_method=source_recon.method;
    S2.dirname=source_recon.dirname;
    S2.use_class_channel=true;
    
    D = osl_inverse_model(S2);
    
    source_recon_results.BF=load([source_recon_sess.dirname '/BF.mat']);

    runcmd(['rm -rf ' source_recon_sess.dirname '/BF.mat']);

    source_recon_results.type = source_recon_sess.type;
    source_recon_results.report=report;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate source recon web report for this session
    
    source_recon_results.report = osl_report_write(source_recon_results.report);        
    source_recon_report         = osl_report_add_sub_report(source_recon_report, source_recon_results.report);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 

    source_recon_results.source_recon=source_recon;
    source_recon_results.mask_fname=mask_fname;
    source_recon_results.mni_coords=mni_coords;        
    source_recon_results.gridstep=source_recon.gridstep;        
    source_recon_results.recon_method=source_recon.method;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save

    source_recon_results.session_name=source_recon.session_name;    
    source_recon_results.fname=[source_recon_results.session_name '_recon' ];
    disp(['Saving source-space results: ' source_recon_results.fname]);  % changed by DM
    
    oat_save_results(oat,source_recon_results);

    results_fnames{sessi}=source_recon_results.fname;
    
    
end;

%%%%%%%%%%%%%%%%%%%
%% summary plots over sessions
source_recon_results.pca_order=nan(length(oat.source_recon.sessions_to_do),1);    
source_recon_results.normalisation=nan(length(source_recon_results.normalisation),length(oat.source_recon.sessions_to_do),1);    
oat.source_recon.results_fnames=results_fnames;
          
for sessi=1:length(oat.source_recon.sessions_to_do), sessnum=oat.source_recon.sessions_to_do(sessi);

    try,
        % load in opt results for this session:            
        res=osl_load_oat_results(oat, oat.source_recon.results_fnames{sessnum});
        
        source_recon_results.pca_order(sessi)=res.pca_order;
        
        for ff=1:length(res.normalisation),
            source_recon_results.normalisation(ff,sessi)=res.normalisation(ff);
        end;

    catch ME,
        disp(['Could not get summary diagnostics for ' oat.source_recon.results_fnames{sessnum}]);
        ME.getReport
    end;
end;

source_recon_report=osl_report_set_figs(source_recon_report,['pca_order']);
plot(oat.source_recon.sessions_to_do,source_recon_results.pca_order,'*');xlabel('sess no.');ylabel(['PCA dim used']); 
source_recon_report=osl_report_print_figs(source_recon_report);

for ff=1:length(res.normalisation),
    source_recon_report=osl_report_set_figs(source_recon_report,[source_recon_results.source_recon.modalities{ff} ' normalisation']);
    plot(oat.source_recon.sessions_to_do,source_recon_results.normalisation(ff,:),'*');xlabel('sess no.');ylabel([source_recon_results.source_recon.modalities{ff} ' normalisation']); 
    source_recon_report=osl_report_print_figs(source_recon_report);
end;

%%%%%%%%%%%%%%%%%%%
%% generate source recon web report
source_recon_report=osl_report_write(source_recon_report);        
source_recon_results.report=source_recon_report;

end
