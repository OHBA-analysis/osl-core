function [source_recon_results]=oat_source_recon_report(oat,source_recon_results)

% report=oat_source_recon_report(oat)
%
% Outputs html page report from OAT source recon stage results
%
% Can adjust behaviour by changing oat.source_recon.report
% settings
%
% MWW 2015

%%%%%%%%%%%%%%%%%%%
%% if needed set diagnostic report up 
if isempty(source_recon_results.report)
    report_dir=[oat.results.plotsdir '/' oat.results.date '_source_recon'];
    source_recon_results.report=osl_report_setup(report_dir,'Source recon');   
end

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
        end
        
        source_recon_results.pca_order(sessi)=res.pca_order;
        
        for ff=1:length(res.normalisation)
            source_recon_results.normalisation(ff,sessi)=res.normalisation(ff);
        end

        %% compute source variance diagnostics
        if oat.source_recon.report.do_source_variance_maps
            D=spm_eeg_load(res.BF.write.spmeeg.files{1});
            [~, mean_V] = osl_source_variance(D);
            source_recon_results.source_variance_nii_fname{sessi}=nii.quicksave(mean_V,[oat.source_recon.dirname '/source_variance_sess' num2str(sessi) '.nii.gz']);

            % diagnostic fig
            S2=[];                
            S2.percrange=[1 99];
            S2.mni_coord=[2 -14 24];
            S2.title='';
            fig_name=['source_variance_sess' num2str(sessi)];  
            fig_title=['Source variance for sess' num2str(sessi) ': fslview(''' source_recon_results.source_variance_nii_fname{sessi} ''')'];
            S2.fname=source_recon_results.source_variance_nii_fname{sessi};
            S2.gridstep=oat.source_recon.gridstep;
            o=ortho_overlay_act( S2 );    
            fig_handle=o.fig; 
            source_recon_results.report=osl_report_set_figs(source_recon_results.report,fig_name,fig_handle,fig_title);        
            source_recon_results.report=osl_report_print_figs(source_recon_results.report);
        end

    catch ME
        disp(['Could not get summary diagnostics for ' oat.source_recon.results_fnames{sessnum}]);
        ME.getReport
    end
end

try
    source_recon_results.report=osl_report_set_figs(source_recon_results.report,'pca_order');
    plot(oat.source_recon.sessions_to_do,source_recon_results.pca_order,'*');xlabel('sess no.');ylabel('PCA dim used'); 
    source_recon_results.report=osl_report_print_figs(source_recon_results.report);
catch
end

try
    for ff=1:length(res.normalisation),
        source_recon_results.report=osl_report_set_figs(source_recon_results.report,[oat.source_recon.modalities{ff} ' normalisation']);
        plot(oat.source_recon.sessions_to_do,source_recon_results.normalisation(ff,:),'*');xlabel('sess no.');ylabel([oat.source_recon.modalities{ff} ' normalisation']); 
        source_recon_results.report=osl_report_print_figs(source_recon_results.report);
    end
catch
end

%%%%%%%%%%%%%%%%%%%
%% generate source recon web report
source_recon_results.report=osl_report_write(source_recon_results.report, oat.results.report);        
disp('View source recon report at:');
disp(['<a href="' source_recon_results.report.html_fname '">' source_recon_results.report.html_fname '</a>']);

end

