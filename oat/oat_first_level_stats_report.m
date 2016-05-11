function report=oat_first_level_stats_report(oat,first_level_results,report_dir_in)

% report=oat_first_level_stats_report(oat,first_level_results)
%
% Outputs html page report from OAT first level stage results
%
% Can adjust behaviour by changing oat.first_level.report
% settings, e.g.
% oat.first_level.report.modality_to_do (only in sensor space)
% oat.first_level.report.first_level_cons_to_do; % plots only these contrasts, 
% and uses first one in list to determine max vox, time, freq
% oat.first_level.report.time_range; % to determine max vox, time, freq
% oat.first_level.report.freq_range; % to determine max vox, time, freq
%
% MWW 2013

global OSLDIR;

if nargin<3,
    report_dir_in='';
end;

% contrasts                                  
contrast_list=oat.first_level.contrast;
for c=1:length(contrast_list),
    contrast_list{c}=contrast_list{c}(:);
end;

if isstr(first_level_results)
    first_level_results=oat_load_results(oat,first_level_results);
else
    first_level_results=first_level_results;
end;

first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts

max_cons=6;
if(length(first_level_cons_to_do)>max_cons)                
    disp(['Only plotting reports for the first ' num2str(max_cons) ' contrasts']);
    first_level_cons_to_do=first_level_cons_to_do(1:max_cons);
end;

times=first_level_results.times;
frequencies=first_level_results.frequencies;
ntpts = length(times);
nfreqs = length(frequencies); 
is_sensor_space=strcmp(first_level_results.recon_method,'none'); % sensor space analysis                    

%%%%%%%%%%%%%%%%%%%
%% set session specific first level diagnostic report up 
if isempty(report_dir_in),
    report_dir=[first_level_results.report.dir '/sess_' first_level_results.session_name];
else
    report_dir=report_dir_in;
end;

report_name=first_level_results.session_name;
if(is_sensor_space)
    report_dir=[report_dir '_' oat.first_level.report.modality_to_do];
    report_name=[report_name ', ' oat.first_level.report.modality_to_do];
end;
report=osl_report_setup(report_dir,report_name);   
   
%%%%%%%%%%%%%%%%%%%
%% diagnostic plot of design matrix

% diagnostic plot of design matrix      
report=osl_report_set_figs(report,'design_matrix',sfigure,'Design Matrix');
title('GLM design matrix');

if(size(first_level_results.cope,2)==1)
    imagesc(first_level_results.D_sensor_data.time(first_level_results.tf_time_indices_into_D_times),1:size(first_level_results.x,2),first_level_results.x');       
    plot4paper('time (secs)','regressor no.');
else
    imagesc(first_level_results.x');
    plot4paper('regressor no.','trial no.');
end;
colorbar;
report=osl_report_print_figs(report);
    
if is_sensor_space,                   

    % find max voxel
    S2=oat.first_level.report; 
    S2.stats=first_level_results;
    S2.modality=oat.first_level.report.modality_to_do;
    [vox_ind_max time_ind_max freq_ind_max stats max_stat] = oat_find_max_stats( S2 );

    if(ntpts>1)
        
        % stats plots at max voxel
        S2=[];
        S2.stats=stats;
        S2.chanlabel=stats.chanlabels{vox_ind_max}; % can change to this setting to plot the time courses at the voxel with the max t-stat
        S2.first_level_cons_to_do=first_level_cons_to_do; % plots all of these contrasts
        [fig_handle fig_name fig_title] = oat_plot_vox_stats(S2);
        fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ')'];
        
        report=osl_report_set_figs(report,fig_name,fig_handle,fig_title);              
        report=osl_report_print_figs(report);
        
    end;
         
    %%%%%%%%%%%%%
    %% need to call oat_stats_multiplotER to get dat and cfg
    clear cfg dats

    S2=[];
    S2.oat=oat;
    S2.stats_fname=first_level_results;
    S2.modality=oat.first_level.report.modality_to_do;
    S2.do_plots=0;

    %%%%%%%%%%%%%
    %% plot topography at time point with max tstat            
    for coni=1:length(first_level_cons_to_do),  
        con=first_level_cons_to_do(coni); 
        S2.first_level_contrast=con;

        fig_handle=[];
        fig_name={};
        fig_title={};
        S2.view_cope=1;
                    
        if(nfreqs>1),
            
            [cfg, dat, fig_handle_tmp]=oat_stats_multiplotTFR(S2);
                    
            cfg.xlim        = [first_level_results.times(time_ind_max) first_level_results.times(time_ind_max)];
            cfg.ylim        = [first_level_results.frequencies(freq_ind_max) first_level_results.frequencies(freq_ind_max)];
            cfg.zlim        = 'maxmin';
            cfg.interactive = 'no';
            cfg.colorbar    = 'yes';   
            
            fig_name{1}=['cope_at_max_multiplot_c' num2str(con)];  
            fig_title{1}=['COPE for c' num2str(con) ' [' first_level_results.first_level_contrast_name{con} ']']; 
            fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ')'];
        
            ft_topoplotTFR(cfg,dat);axis tight;
            title([oat.first_level.report.modality_to_do]);
            
            report=osl_report_set_figs(report,fig_name,fig_handle,fig_title);                    
            report=osl_report_print_figs(report);
        else
            S2.do_plots=1;
            
            fig_name{1}=['cope_at_max_multiplot_c' num2str(con)];  
            fig_title{1}=['COPE for c' num2str(con) ' [' first_level_results.first_level_contrast_name{con} ']']; 
            fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ')'];
        
            [cfg, dat, fig_handle(1)]=oat_stats_multiplotER(S2);
            view([90 90])
            
            fig_handle(2)=sfigure; 
            fig_name{2}=['cope_at_max_topo_c' num2str(con)];  
            fig_title{2}=['COPE for c' num2str(con) ' [' first_level_results.first_level_contrast_name{con} ']']; 
            fig_title{2}=[fig_title{2} ' (Max found using c' num2str(first_level_cons_to_do(1)) ')'];
        
            cfg.xlim        = [first_level_results.times(time_ind_max) first_level_results.times(time_ind_max)];
            cfg.ylim        = 'maxmin';
            cfg.zlim        = 'maxmin';
            cfg.interactive = 'no';
            cfg.colorbar    = 'yes';
            
            ft_topoplotTFR(cfg,dat{1});axis tight;
            title([oat.first_level.report.modality_to_do]);
            
            report=osl_report_set_figs(report,fig_name,fig_handle,fig_title);                    
            report=osl_report_print_figs(report);
        end;
              
    end;                           

else
    
    % source space
    
    % find max voxel
    S2=oat.first_level.report;    
    S2.stats=first_level_results;               
    [vox_ind_max time_ind_max freq_ind_max stats] = oat_find_max_stats( S2 );

    if(ntpts>1)       
        % stats plots at max voxel
        S2=[];
        S2.stats=stats;
        S2.vox_coord=stats.mni_coords(vox_ind_max,:); % can change to this setting to plot the time courses at the voxel with the max t-stat
        S2.first_level_cons_to_do=first_level_cons_to_do; % plots all of these contrasts
        [fig_handle fig_name fig_title] = oat_plot_vox_stats(S2);   
        fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ')'];        
        report=osl_report_set_figs(report,fig_name,fig_handle,fig_title);        
        report=osl_report_print_figs(report);
    end;
    
    if ~isfield(oat.first_level, 'mni_coords'),

        %%%%%%%%%%%%%%%%%%%
        %% output required niis
        S2=[];
        S2.oat=oat;
        if(ntpts>1)
            S2.time_range=[first_level_results.times(time_ind_max) first_level_results.times(time_ind_max)];
        end;
        S2.stats=first_level_results;
        %S2.stats_fname=oat.first_level.results_fnames{1};
        S2.first_level_contrasts=first_level_cons_to_do; % list of first level contrasts to output
        resamp_gridstep=2;
        S2.freq_bin=freq_ind_max;
        S2.resamp_gridstep=resamp_gridstep;%oat.source_recon.gridstep;
        [statsdir,times,count]=oat_save_nii_stats(S2);

        %%%%%%%%%%%%%
        %% plot cope orthoviews at time point (and freq) with max tstat 
        fig_handle=[];
        fig_name={};
        fig_title={}; 
        S2=[];                
        
        S2.percrange=[96 99.9];
        S2.mni_coord=first_level_results.mni_coords(vox_ind_max,:);
        S2.title='';

        for coni=1:length(first_level_cons_to_do),  
            con=first_level_cons_to_do(coni); 

            fig_handle(coni)=sfigure; 
            fig_name{coni}=['cope_at_maxt_smap_c' num2str(con)];  
            fig_title{coni}=['COPE for c' num2str(con) ' [' first_level_results.first_level_contrast_name{con} '] at t=' num2str(first_level_results.times(time_ind_max)) 'secs' ', f=' mat2str(first_level_results.frequency_ranges(freq_ind_max,:),3) 'Hz'];                 
            fig_title{coni}=[fig_title{coni} ' (Max found using c' num2str(first_level_cons_to_do(1)) ')'];
            set(fig_handle(coni),'Position',[1 1 1300 450]);
            S2.fname=[statsdir '/cope' num2str(con) '_' num2str(resamp_gridstep) 'mm'];                              
            ortho_overlay_act( S2 );                     

        end;
    
    end;

    report=osl_report_set_figs(report,fig_name,fig_handle,fig_title);        
    report=osl_report_print_figs(report);

end;    

%%%%%%%%%%%%%%%%%%%
%% generate source recon web report
report=osl_report_write(report);        
disp(['View first level report at:']);
disp([report.html_fname]);

end

