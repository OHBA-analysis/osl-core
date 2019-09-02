function [results_fnames source_recon_results]=osl_run_source_recon_sensorspace(oat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup sensor space data
%
% results_fnames=osl_run_source_recon_sensorspace(oat)
%
% Mark Woolrich 2012

getenv('OSLDIR');

source_recon=oat.source_recon;

source_recon_results=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set first level diagnostic report up    
report_dir=[oat.results.plotsdir '/' oat.results.date '_source_recon'];
source_recon_results.report=osl_report_setup(report_dir,['Source recon - sensor space data setup']);   

for sessi_todo=1:length(oat.source_recon.sessions_to_do),   

    sessi=oat.source_recon.sessions_to_do(sessi_todo);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OAT SOURCE RECON (SENSOR SPACE SETUP) ON SESS = ' num2str(sessi) '  %%%%%%%%%%%%%%%%%%%%%%%'])
    
    % setup session specific diagnostic report up
    report_dir = [source_recon_results.report.dir '/sess' num2str(sessi)];
    report_sess     = osl_report_setup(report_dir,['Session ' num2str(sessi)]); 
    
    source_recon_sess=source_recon;
    
    source_recon_sess.do_plots=oat.do_plots;
    source_recon_sess.session_name=['session' num2str(sessi)];

    if ~isempty(source_recon.D_continuous),        
        [p fname e] = fileparts(source_recon.D_continuous{sessi});       
        source_recon_sess.D_continuous=[p '/' fname '.mat'];       
        disp('Using continuous data as input');        
    else
        source_recon_sess.D_continuous=[];
    end;
    
    if ~isempty(source_recon.D_epoched),
        [p fname e] = fileparts(source_recon.D_epoched{sessi});       
        source_recon_sess.D_epoched=[p '/' fname '.mat'];       
        disp('Using epoched data as input');
    else
        source_recon_sess.D_epoched=[];
    end;    
    
%    source_recon_sess.mri=source_recon.mri{sessi};
    
    if length(source_recon.pca_dim)>1,
        source_recon_sess.pca_dim=source_recon.pca_dim(sessi);
    else
        source_recon_sess.pca_dim=source_recon.pca_dim;        
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    [source_recon_sess source_recon_sess_results source_recon_sess_results.report] = oat_prepare_source_recon(source_recon_sess, report_sess);
    
    %%%%%%%%%%%%%%%%%%%
    %% generate source recon web report for this session
    source_recon_sess_results.report = osl_report_write(source_recon_sess_results.report, source_recon_results.report);         
    source_recon_results.report         = osl_report_add_sub_report(source_recon_results.report, source_recon_sess_results.report);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     
    
    source_recon_sess_results.source_recon=source_recon;
    source_recon_sess_results.recon_method=source_recon.method;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save
    
    source_recon_sess_results.session_name=source_recon_sess.session_name;    
    source_recon_sess_results.fname=[source_recon_sess_results.session_name '_recon' ];
    disp(['Saving source recon (sensorspace setup) results: ' source_recon_sess_results.fname]);
    
    oat_save_results(oat,source_recon_sess_results);

    results_fnames{sessi}=source_recon_sess_results.fname;    

end

%%%%%%%%%%%%%%%%%%%
%% summary plots over sessions
source_recon_results.pca_order=nan(length(oat.source_recon.sessions_to_do),1);    
source_recon_results.normalisation=nan(length(oat.source_recon.sessions_to_do),1);    
        
oat.source_recon.results_fnames=results_fnames;

for sessi=1:length(oat.source_recon.sessions_to_do), sessnum=oat.source_recon.sessions_to_do(sessi);

    try
        % load in oat results for this session:            
        res=oat_load_results(oat, oat.source_recon.results_fnames{sessnum});

        mod_ind=find(strcmp(oat.source_recon.modalities,oat.first_level.report.modality_to_do));                
        
        for ff=1:length(res.pca_order)
            source_recon_results.pca_order(ff,sessi)=res.pca_order(ff);
        end
        for ff=1:length(res.normalisation)
            source_recon_results.normalisation(ff,sessi)=res.normalisation(ff);
        end

    catch ME
        disp(['Could not get summary diagnostics for ' oat.source_recon.results_fnames{sessnum}]);
        ME.getReport
    end
end

for ff=1:length(res.pca_order)
    source_recon_results.report=osl_report_set_figs(source_recon_results.report,[oat.source_recon.modalities{ff} ' pca_order']);
    plot(oat.source_recon.sessions_to_do,source_recon_results.pca_order(ff,:),'*');xlabel('sess no.');ylabel([oat.source_recon.modalities{ff} ' PCA dim']); 
    source_recon_results.report=osl_report_print_figs(source_recon_results.report);
end

for ff=1:length(res.normalisation)
    source_recon_results.report=osl_report_set_figs(source_recon_results.report,[oat.source_recon.modalities{ff} ' normalisation']);
    plot(oat.source_recon.sessions_to_do,source_recon_results.normalisation(ff,:),'*');xlabel('sess no.');ylabel([oat.source_recon.modalities{ff} ' normalisation']); 
    source_recon_results.report=osl_report_print_figs(source_recon_results.report);
end

%%%%%%%%%%%%%%%%%%%
%% generate source recon web report
source_recon_results.report=osl_report_write(source_recon_results.report, oat.results.report);        

end
