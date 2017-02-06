function [report, source_recon_results]=oat_source_recon_report(oat,report,source_recon_results)

% report=oat_source_recon_report(oat)
%
% Outputs html page report from OAT source recon stage results
%
% Can adjust behaviour by changing oat.source_recon.report
% settings
%
% MWW 2015

if nargin<2,
    report=[];
end;

if nargin<3,
    source_recon_results=[];
end;

%%%%%%%%%%%%%%%%%%%
%% if needed set diagnostic report up 
if isempty(report),
    report_dir=[oat.results.plotsdir '/' oat.results.date '_source_recon'];
    report=osl_report_setup(report_dir,'Source recon');   
end;

%%%%%%%%%%%%%%%%%%%
%% setup containers
source_recon_results.pca_order=nan(length(oat.source_recon.sessions_to_do),1);    
  
%%%%%%%%%%%%%%%%%%%
%% loop over sessions 
for sessi=1:length(oat.source_recon.sessions_to_do), sessnum=oat.source_recon.sessions_to_do(sessi);

    try
        %% load in opt results for this session:            
        res=oat_load_results(oat, oat.source_recon.results_fnames{sessnum});
        
        if sessi==1
            source_recon_results.normalisation=nan(length(res.normalisation),length(oat.source_recon.sessions_to_do),1);    
        end;
        
        source_recon_results.pca_order(sessi)=res.pca_order;
        
        for ff=1:length(res.normalisation),
            source_recon_results.normalisation(ff,sessi)=res.normalisation(ff);
        end;

        %% compute source variance diagnostics
        if oat.source_recon.report.do_source_variance_maps
            D=spm_eeg_load(res.BF.write.spmeeg.files{1});
            [~, mean_V] = osl_source_variance(D);
            source_recon_results.source_variance_nii_fname{sessi}=nii_quicksave(mean_V,[oat.source_recon.dirname '/source_variance_sess' num2str(sessi)]);

            % diagnostic fig
            S2=[];                
            S2.percrange=[1 99];
            S2.mni_coord=[2 -14 24];
            S2.title='';
            fig_handle=sfigure; 
            fig_name=['source_variance_sess' num2str(sessi)];  
            fig_title=['Source variance for sess' num2str(sessi) ': fslview(''' source_recon_results.source_variance_nii_fname{sessi} ''')'];
            set(fig_handle,'Position',[1 1 1300 450]);
            S2.fname=source_recon_results.source_variance_nii_fname{sessi};
            S2.gridstep=oat.source_recon.gridstep;
            ortho_overlay_act( S2 );    
            report=osl_report_set_figs(report,fig_name,fig_handle,fig_title);        
            report=osl_report_print_figs(report);
        end;

    catch ME,
        disp(['Could not get summary diagnostics for ' oat.source_recon.results_fnames{sessnum}]);
        ME.getReport
    end
end

try
    report=osl_report_set_figs(report,'pca_order');
    plot(oat.source_recon.sessions_to_do,source_recon_results.pca_order,'*');xlabel('sess no.');ylabel('PCA dim used'); 
    report=osl_report_print_figs(report);
catch
end

try
    for ff=1:length(res.normalisation),
        report=osl_report_set_figs(report,[oat.source_recon.modalities{ff} ' normalisation']);
        plot(oat.source_recon.sessions_to_do,source_recon_results.normalisation(ff,:),'*');xlabel('sess no.');ylabel([oat.source_recon.modalities{ff} ' normalisation']); 
        report=osl_report_print_figs(report);
    end
catch
end

%%%%%%%%%%%%%%%%%%%
%% generate source recon web report
report=osl_report_write(report);        
disp('View source recon report at:');
disp(['<a href="' report.html_fname '">' report.html_fname '</a>']);

end

