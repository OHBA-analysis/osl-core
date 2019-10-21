function [results_fnames, source_recon_results]=oat_run_source_recon(oat)

% results_fnames=osl_run_source_recon_beamform(oat)
%
% takes in an OAT, which needs to be setup by calling oat=osl_setup_oat(S), struct 
% and runs beamformer 
% 
% This function should normally be called using osl_run_oat(oat);
%
% MW 2011

OSLDIR = getenv('OSLDIR');

source_recon=oat.source_recon;

source_recon_results=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup mni coords using std space mask
if(isfield(source_recon,'mask_fname'))

    if ~isempty(source_recon.gridstep)
        mask_fname=nii.resample(source_recon.mask_fname, ...
            [source_recon.mask_fname '_' num2str(source_recon.gridstep) 'mm.nii.gz'], ...
            source_recon.gridstep, 'interptype', 'nearest','enforce_mask',true);
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
    mask_fname=[osldir '/std_masks/MNI152_T1_' num2str(source_recon.gridstep) 'mm_brain.nii.gz']
    disp(['Using whole brain']);
    
    % setup mni coords usource_recong std space mask
    [ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname);    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save mask in oat dir
if(~strcmp(mask_fname,''))
   copyfile(mask_fname,fullfile(oat.source_recon.dirname,'source_recon_mask.nii.gz'));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set diagnostic report up    
report_dir=[oat.results.plotsdir '/' oat.results.date '_source_recon'];
source_recon_results.report=osl_report_setup(report_dir,['Source recon']);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp(['Using ' source_recon.method ' for source reconstruction']);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over sessions
for sessi_todo=1:length(source_recon.sessions_to_do),   
        
    sessi = source_recon.sessions_to_do(sessi_todo);
 
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OAT SOURCE RECON ON SESS = ' num2str(sessi) '  %%%%%%%%%%%%%%%%%%%%%%%'])
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set session specific diagnostic report up    
    report_dir = [source_recon_results.report.dir '/sess' num2str(sessi)];
    report_sess     = osl_report_setup(report_dir,['Session ' num2str(sessi)]);       

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
        
    if length(source_recon.pca_dim)>1,
        source_recon_sess.pca_dim=source_recon.pca_dim(sessi);
    else
        source_recon_sess.pca_dim=source_recon.pca_dim;        
    end;
    
    source_recon_sess.mni_coords = mni_coords;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Prepare data (bandpass filter, epoch, normalise sensors)
    [source_recon_sess source_recon_sess_results source_recon_sess_results.report] = oat_prepare_source_recon(source_recon_sess, report_sess);
    D=source_recon_sess_results.BF.data.D;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DO CO-REGISTRATION AND HEAD MODEL
    % Need to re-do it now following any montaging done in osl_prepare_oat_source_recon

    if ~strcmp(oat.source_recon.normalise_method,'none')
        
        disp('Co-registration and setting up forward model...');

        D = osl_forward_model(D,source_recon_sess);
        D.save();

        if(oat.do_plots)
            %spm_eeg_inv_checkmeshes(D);
            spm_eeg_inv_checkdatareg(D);
            %spm_eeg_inv_checkforward(D, 1);
        end;

    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run inverse
    
    S2=[];
    S2.timespan=source_recon_sess.time_range;
    S2.type=source_recon_sess.type;
    S2.pca_order=source_recon_sess.pca_dim;
    S2.modalities=source_recon_sess.modalities;
    S2.fuse='meg';
    S2.conditions=source_recon_sess.conditions;
    S2.inverse_method=source_recon_sess.method;
    S2.dirname=source_recon_sess.dirname;
    S2.dirname = D.path;
    S2.use_class_channel=true;
    D = osl_inverse_model(D,mni_coords,S2);
    
    source_recon_sess_results.BF=load(fullfile(source_recon_sess.dirname,'BF.mat'));

    delete(fullfile(source_recon_sess.dirname,'BF.mat'));

    source_recon_sess_results.type = source_recon_sess.type;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate source recon web report for this session
    
    source_recon_sess_results.report = osl_report_write(source_recon_sess_results.report, source_recon_results.report);         
    source_recon_results.report         = osl_report_add_sub_report(source_recon_results.report, source_recon_sess_results.report);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 

    source_recon_sess_results.source_recon=source_recon;
    source_recon_sess_results.mask_fname=mask_fname;
    source_recon_sess_results.mni_coords=mni_coords;        
    source_recon_sess_results.gridstep=source_recon.gridstep;        
    source_recon_sess_results.recon_method=source_recon.method;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save

    source_recon_sess_results.session_name=source_recon.session_name;    
    source_recon_sess_results.fname=[source_recon_sess_results.session_name '_recon' ];
    disp(['Saving source-space results: ' source_recon_sess_results.fname]);  
    
    oat_save_results(oat,source_recon_sess_results);

    results_fnames{sessi}=source_recon_sess_results.fname;
    
end

oat.source_recon.results_fnames=results_fnames;

%%%%%%%%%%%%%%%%%%%
%% summary plots over sessions
source_recon_results = oat_source_recon_report(oat, source_recon_results);

end
