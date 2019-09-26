function group_level_report = oat_group_level_stats_report(oat,group_level_results,report_dir_in)

% report = oat_group_level_stats_report(oat,group_level_results)
%
% Outputs html page report from OAT group stage results
%
% Adjust behaviour by changing oat.group_level.report
% settings, e.g.
% oat.group_level.report.modality_to_do (only in sensor space)
% oat.group_level.report.group_level_cons_to_do; % plots all of these contrasts, and uses first one in list to determine max vox, time, freq
% oat.group_level.report.first_level_cons_to_do; % plots all of these contrasts, and uses first one in list to determine max vox, time, freq
% oat.group_level.report.time_range; % to determine max vox, time, freq
% oat.group_level.report.freq_range; % to determine max vox, time, freq
% oat.group_level.report.show_lower_level_cope_maps
% oat.group_level.report.show_lower_level_copes
%
% MWW 2013

OSLDIR = getenv('OSLDIR');

if nargin<3
    report_dir_in='';
end

if isstr(group_level_results)
    group_level_results=oat_load_results(oat,group_level_results);
else
    group_level_results=group_level_results;
end

if(~isempty(setdiff(oat.group_level.report.first_level_cons_to_do,oat.group_level.first_level_contrasts_to_do))),
    error('oat.group_level.report.first_level_cons_to_do must be a subset of oat.group_level.first_level_contrasts_to_do');
end

if(max(oat.group_level.report.group_level_cons_to_do)>length(oat.group_level.group_contrast)),
    error('oat.group_level.report.group_cons_to_do must be a subset of oat.group_level.group_contrasts');
end

group_level_cons_to_do=oat.group_level.report.group_level_cons_to_do; % plots all of these contrasts, and uses first one in list to determine max vox, time, freq
first_level_cons_to_do=oat.group_level.report.first_level_cons_to_do; % plots all of these contrasts, and uses first one in list to determine max vox, time, freq

max_cons=6;
if(length(first_level_cons_to_do)>max_cons)                
    disp(['Only plotting reports for the first ' num2str(max_cons) ' contrasts']);
    first_level_cons_to_do=first_level_cons_to_do(1:max_cons);
end

times=group_level_results.times;
frequencies=group_level_results.frequencies;
ntpts = length(times);
nfreqs = length(frequencies); 
is_sensor_space=strcmp(oat.source_recon.method,'none'); % sensor space analysis                    
    
%%%%%%%%%%%%%%%%%%%
%% set first level diagnostic report up    
if isempty(report_dir_in),
    report_dir=[oat.results.plotsdir '/' oat.first_level.name '_' oat.subject_level.name '_' oat.group_level.name];
else
    report_dir=report_dir_in;
end

report_name=['Group Level'];
if(is_sensor_space)
    report_dir=[report_dir '_' oat.group_level.report.modality_to_do];
    report_name=[report_name ', ' oat.group_level.report.modality_to_do];
end

group_level_report=osl_report_setup(report_dir,report_name);   

%%%%%%%%%%%%%%%%%%%
%% diagnostic plot of group design matrix

group_level_report=osl_report_set_figs(group_level_report,'design_matrix',sfigure,'Design Matrix');
imagesc(oat.group_level.group_design_matrix');
title('GLM design matrix');plot4paper('regressor no.','subject no.');
group_level_report=osl_report_print_figs(group_level_report);


% find max voxel, time, freq for first first level con and first group
% level con in the list of contrasts to show
S2=oat.group_level.report;
S2.stats=group_level_results;
S2.modality=oat.group_level.report.modality_to_do; % in case we are in sensor space
[vox_ind_max time_ind_max freq_ind_max stats max_stat] = oat_find_max_stats( S2 );

%%%%%%%%%%%%%%%%%%%
%% lower level cope plots
if(oat.group_level.report.show_lower_level_copes)
     if 0,
        fig_handle=[];
        fig_name={};
        fig_title={};
        gcon = oat.group_level.report.group_level_cons_to_do(1); % only show for first group  level con in the list of contrasts to show
        for coni=1:length(first_level_cons_to_do),  
            con=first_level_cons_to_do(coni); 

            fig_handle(coni)=sfigure; 
            set(fig_handle(coni),'Position',[1 1 1300 450]);          
            fig_name{coni}=['lower_level_stdcope' num2str(con) '_gc' num2str(gcon)];  
            fig_title{coni}=['Average (over time, freq and space) Lower level STDCOPEs for c' num2str(con) ' [' oat.first_level.contrast_name{con} ']' ', gc' num2str(gcon)];  

            % lower_level_stdcopes is nchans x nsubjs x ntpts x nfreqs
            plot(permute(mean(mean(mean(group_level_results.lower_level_stdcopes{con},1),3),4),[2,1,3,4]),'LineWidth',2);ho;
            plot4paper('subject','1st level std cope'); 

        end
        group_level_report=osl_report_set_figs(group_level_report,fig_name,fig_handle,fig_title);        
        group_level_report=osl_report_print_figs(group_level_report);
     end
end


%%%%%%%%%%%%%%%%%%%
%% sub report for each first level contrast    
for cc=1:length(first_level_cons_to_do),

    con4report=first_level_cons_to_do(cc);

    disp(['Processing Group stats for first level contrast ' num2str(con4report) ' [' oat.first_level.contrast_name{con4report} ']']);
    disp(['Number ' num2str(cc) ' of ' num2str(length(oat.group_level.first_level_contrasts_to_do))]);

    %%%%%%%%%%%%%%%%%%%
    %% set first level con diagnostic report up    
    con_report_dir=[group_level_report.dir '/c' num2str(con4report) ];
    con_report=osl_report_setup(con_report_dir,['First level contrast ' num2str(con4report) ' [' oat.first_level.contrast_name{con4report} ']']); 

    if(ntpts>1)

        % find max voxel, time, freq
        S2=oat.group_level.report;
        S2.stats=group_level_results;
        first_level_cons_to_do_reordered_for_con=[con4report, setdiff(first_level_cons_to_do,first_level_cons_to_do(con4report))];
        S2.first_level_cons_to_do=first_level_cons_to_do_reordered_for_con;
        S2.modality=oat.group_level.report.modality_to_do; % in case we are in sensor space

        [vox_ind_max time_ind_max freq_ind_max stats max_stat] = oat_find_max_stats( S2 );

        if ~is_sensor_space
            spatial_loc_str=mat2str(stats.mni_coords(vox_ind_max,:),4);            
        else
            spatial_loc_str=[stats.chanlabels{vox_ind_max}];
        end
        
        % stats plots at max voxel
        % show for ALL group contrasts
        
        for gcc=1:length(group_level_cons_to_do)

            gccon=group_level_cons_to_do(gcc);

            S2=[];
            S2.stats=stats;
            if ~is_sensor_space
                S2.vox_coord=stats.mni_coords(vox_ind_max,:);
            else
                S2.chanlabel=stats.chanlabels{vox_ind_max};
            end
            S2.first_level_cons_to_do=first_level_cons_to_do; % plots all of these contrasts
            S2.group_level_cons_to_do=gccon; % plots all of these contrasts
            [fig_handle fig_name fig_title] = oat_plot_vox_stats(S2);
            fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do_reordered_for_con(1)) ', gc' num2str(gccon) ')'];
            con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
            con_report=osl_report_print_figs(con_report);
        end
        
        %%%%%%%%%%%%%%%%%%%
        %% lower level cope plots: these will be the same regardless of the group con
        if(oat.group_level.report.show_lower_level_copes && length(group_level_results.frequencies)==1)
            if(isempty(oat.group_level.time_range))    
                oat.group_level.time_range=[group_level_results.lower_level_times(1) group_level_results.lower_level_times(end)];
            end

            fig_handle=[];
            fig_name={};
            fig_title={};

            time_idx = group_level_results.lower_level_times > oat.group_level.time_range(1) & ...
                         group_level_results.lower_level_times < oat.group_level.time_range(2);
            time_range = group_level_results.lower_level_times(time_idx);

            for coni=1:length(first_level_cons_to_do)  
                con=first_level_cons_to_do(coni);
                fig_handle(coni)=sfigure; 
                set(fig_handle(coni),'Position',[1 1 1300 450]);          
                fig_name{coni}=['lower_level_cope_ts_at_maxt_c' num2str(con) '_gc' num2str(con)];  
                fig_title{coni}=['Lower level COPEs for c' num2str(con) ' [' oat.first_level.contrast_name{con} '], at ' spatial_loc_str];  
                fig_title{coni}=[fig_title{coni} ' (Max found using c' num2str(first_level_cons_to_do_reordered_for_con(1)) ', gc' num2str(gccon) ')'];

                plot(time_range,permute(group_level_results.lower_level_copes{con}(vox_ind_max,:,time_idx),[2,3,1,4]),'LineWidth',2);ho;
                plot(time_range,squeeze(mean(group_level_results.lower_level_copes{con}(vox_ind_max,:,time_idx),2)),'LineWidth',6);
                plot4paper('time (secs)','1st level cope'); 
            end

            con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
            con_report=osl_report_print_figs(con_report);

            %% lower level tstats time series
            fig_handle=[];
            fig_name={};
            fig_title={};

            for coni=1:length(first_level_cons_to_do),    
                con=first_level_cons_to_do(coni);
                fig_handle(coni)=sfigure; 
                set(fig_handle(coni),'Position',[1 1 1300 450]);          
                fig_name{coni}=['lower_level_tstats_ts_at_maxt_c' num2str(con) '_gc' num2str(con)];  
                fig_title{coni}=['Lower level tstats for c' num2str(con) ' [' oat.first_level.contrast_name{con} '], at ' spatial_loc_str];
                fig_title{coni}=[fig_title{coni} ' (Max found using c' num2str(first_level_cons_to_do_reordered_for_con(1)) ', gc' num2str(gccon) ')'];
                %subplot(1,2,1);
                plot(group_level_results.times,permute(group_level_results.lower_level_copes{con}(vox_ind_max,:,:)./group_level_results.lower_level_stdcopes{con}(vox_ind_max,:,:),[2,3,1]),'LineWidth',2);ho;
                plot(group_level_results.times,squeeze(mean(group_level_results.lower_level_copes{con}(vox_ind_max,:,:)./group_level_results.lower_level_stdcopes{con}(vox_ind_max,:,:),2)),'LineWidth',6);
                plot4paper('time (secs)','1st level tstats'); 
            end

            con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
            con_report=osl_report_print_figs(con_report);

            if false
                %% lower level std cope
                for coni=1:length(first_level_cons_to_do)   
                    con=first_level_cons_to_do(coni);
                    fig_handle(coni)=sfigure; 
                    set(fig_handle(coni),'Position',[1 1 1300 450]);          
                    fig_name{coni}=['lower_level_stdcope_ts_at_maxt_c' num2str(con) '_gc' num2str(con)];  
                    fig_title{coni}=['Lower level STDCOPEs for c' num2str(con) ' [' oat.first_level.contrast_name{con} '], at ' spatial_loc_str];
                    fig_title{coni}=[fig_title{coni} ' (Max found using c' num2str(first_level_cons_to_do_reordered_for_con(1)) ', gc' num2str(gccon) ')'];
                    %subplot(1,2,1);
                    plot(permute(mean(group_level_results.lower_level_stdcopes{con}(vox_ind_max,:,:),3),[2,3,1]),'LineWidth',2);ho;
                    plot4paper('subject','1st level std cope'); 
                end
                con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
                con_report=osl_report_print_figs(con_report);
            end

        end
        
    end
        
    if is_sensor_space                 
   
        % sensor space

        if ~oat.group_level.space_average,

            %%%%%%%%%%%%%
            %% plot topography at time point with max tstat            

            fig_handle=[];
            fig_name={};
            fig_title={};

            for gconi=1:length(group_level_cons_to_do),  
                gcon=group_level_cons_to_do(gconi); 

                % DO COPEs
                S2=[];
                S2.oat=oat;
                S2.stats_fname=group_level_results;
                S2.modality=oat.group_level.report.modality_to_do;
                S2.do_plots=0;

                S2.first_level_contrast=con;                        
                S2.group_level_contrast=gcon;

                if(nfreqs>1),                
                    % need to call oat_stats_multiplot to get dat and cfg
                    S2.view_cope=1;
                    [cfg, dat, fig_handle_tmp]=oat_stats_multiplotTFR(S2);      
                    S2.view_cope=0;
                    [cfg_tstat, dat_tstat, fig_handle_tmp]=oat_stats_multiplotTFR(S2);      

                    fig_handle(1)=sfigure; 
                    fig_name{1}=['stats_at_max_topo_c' num2str(con4report) '_gc' num2str(gcon)];                      
                    fig_title{1}=['Stats for gc' num2str(gcon) ' [' oat.group_level.group_contrast_name{gcon} ']']; 
                    fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ', gc' num2str(group_level_cons_to_do(1)) ')'];

                    snugplot(1,2,1);
                    cfg.xlim        = [group_level_results.times(time_ind_max) group_level_results.times(time_ind_max)];
                    cfg.ylim        = [group_level_results.frequencies(freq_ind_max) group_level_results.frequencies(freq_ind_max)];
                    cfg.zlim        = 'maxmin';
                    cfg.interactive = 'no';
                    cfg.colorbar    = 'yes';

                    ft_topoplotTFR(cfg,dat);axis tight;
                    title(['COPE']);                    

                    snugplot(1,2,2);
                    cfg_tstat.xlim        = [group_level_results.times(time_ind_max) group_level_results.times(time_ind_max)];
                    cfg_tstat.ylim        = [group_level_results.frequencies(freq_ind_max) group_level_results.frequencies(freq_ind_max)];
                    cfg_tstat.zlim        = 'maxmin';
                    cfg_tstat.interactive = 'no';
                    cfg_tstat.colorbar    = 'yes';

                    ft_topoplotTFR(cfg_tstat,dat_tstat);axis tight;
                    title(['COPE']);

                    con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
                    con_report=osl_report_print_figs(con_report);

                else,

                    
                    
                    S2.do_plots=1;
                    S2.view_cope=1;
                    % need to call oat_stats_multiplot to get dat and cfg                        
                    [cfg, dat, fig_handle(1)]=oat_stats_multiplotER(S2);
                    view([90 90])
                    fig_name{1}=['cope_at_max_multiplot_c' num2str(con4report) '_gc' num2str(gcon)];                      
                    fig_title{1}=['COPE for gc' num2str(gcon) ' [' oat.group_level.group_contrast_name{gcon} '], at t=' num2str(group_level_results.times(time_ind_max)) 'secs, f=' mat2str(group_level_results.frequency_ranges(freq_ind_max,:)) 'Hz'];  
                    fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ', gc' num2str(group_level_cons_to_do(1)) ')'];
                    
                    con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
                    con_report=osl_report_print_figs(con_report);

                    S2.do_plots=0;
                    S2.view_cope=0;
                    [cfg_tstat, dat_tstat, fig_handle_tmp]=oat_stats_multiplotER(S2);

                    fig_handle(1)=sfigure;                    
                    fig_name{1}=['tstats_at_max_topo_c' num2str(con4report) '_gc' num2str(gcon)];                      
                    fig_title{1}=['T-stat for gc' num2str(gcon) ' [' oat.group_level.group_contrast_name{gcon} '], at t=' num2str(group_level_results.times(time_ind_max)) 'secs, f=' mat2str(group_level_results.frequency_ranges(freq_ind_max,:)) 'Hz'];  
                    fig_title{1}=[fig_title{1} ' (Max found using c' num2str(first_level_cons_to_do(1)) ', gc' num2str(group_level_cons_to_do(1)) ')'];
                    
                    snugplot(1,2,1);
                    cfg.xlim        = [group_level_results.times(time_ind_max) group_level_results.times(time_ind_max)];
                    cfg.ylim        = 'maxmin';
                    cfg.zlim        = 'maxmin';
                    cfg.interactive = 'no';
                    cfg.colorbar    = 'yes';                    
                    ft_topoplotTFR(cfg,dat{1});axis tight;
                    title(['COPE']);

                    snugplot(1,2,2);
                    cfg_tstat.xlim        = [group_level_results.times(time_ind_max) group_level_results.times(time_ind_max)];
                    cfg_tstat.ylim        = 'maxmin';
                    cfg_tstat.zlim        = 'maxmin';
                    cfg_tstat.interactive = 'no';
                    cfg_tstat.colorbar    = 'yes';                    
                    ft_topoplotTFR(cfg_tstat,dat_tstat{1});axis tight;
                    title(['T-stat']);

                    con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
                    con_report=osl_report_print_figs(con_report);

                end

            end                         

        end

    else
        
        % source space
        
        if ~oat.group_level.space_average

            %%%%%%%%%%%%%%%%%%%
            %% output required niis
            S2=[];
            S2.oat=oat;
            S2.time_range=[group_level_results.times(time_ind_max) group_level_results.times(time_ind_max)];
            S2.stats=group_level_results;
            %S2.stats_fname=oat.first_level.results_fnames{1};
            S2.first_level_contrasts=[con4report]; % list of first level contrasts to output
            S2.group_level_contrasts = oat.group_level.report.group_level_cons_to_do;
            resamp_gridstep=2;
            S2.resamp_gridstep=resamp_gridstep;%oat.source_recon.gridstep;
            S2.freq_bin=freq_ind_max;  
            S2.stats_dir=report_dir;
  
            [statsdir,times,count]=oat_save_nii_stats(S2);

            %%%%%%%%%%%%%
            %% plot cope orthoviews at time point (and freq) with max tstat 
            fig_handle=[];
            fig_name={};
            fig_title={};
            S2=[];                
            S2.percrange=[96 99.9];
            %S2.percrange=[70 99.9];
            S2.mni_coord=group_level_results.mni_coords(vox_ind_max,:);
            S2.title='';

            for gconi=1:length(group_level_cons_to_do),  
                gcon=group_level_cons_to_do(gconi); 
                fig_name{gconi}=['cope_at_maxt_smap_c' num2str(con4report) '_gc' num2str(gcon)];  
                fig_title{gconi}=['COPE for gc' num2str(gcon) ' [' oat.group_level.group_contrast_name{gcon} '], at t=' num2str(group_level_results.times(time_ind_max)) 'secs, f=' mat2str(group_level_results.frequency_ranges(freq_ind_max,:)) 'Hz'];  
                fig_title{gconi}=[fig_title{gconi} ' (Max found using c' num2str(first_level_cons_to_do(1)) ', gc' num2str(group_level_cons_to_do(1)) ')'];
                S2.fname=[statsdir '/cope' num2str(con4report) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'];                              
                o = ortho_overlay_act( S2 );
                fig_handle(gconi)= o.fig; 
                setpixelposition(fig_handle(gconi),[1 1 1300 450]);
            end
            con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
            con_report=osl_report_print_figs(con_report);
            
            
            %%%%%%%%%%%%%
            %% plot tstat orthoviews at time point with max tstat            
            fig_handle=[];
            fig_name={};
            fig_title={};
            for gconi=1:length(group_level_cons_to_do),  
                gcon=group_level_cons_to_do(gconi); 

                fig_name{gconi}=['tstat_at_maxt_smap_c' num2str(con4report) '_gc' num2str(gcon)];  
                fig_title{gconi}=['T-stat for gc' num2str(gcon) ' [' oat.group_level.group_contrast_name{gcon} '], at t=' num2str(group_level_results.times(time_ind_max)) 'secs, f=' mat2str(group_level_results.frequency_ranges(freq_ind_max,:)) 'Hz'];  
                fig_title{gconi}=[fig_title{gconi} ' (Max found using c' num2str(first_level_cons_to_do(1)) ', gc' num2str(group_level_cons_to_do(1)) ')'];
                S2.fname=[statsdir '/tstat' num2str(con4report) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'];                              
                o = ortho_overlay_act( S2 );
                fig_handle(gconi)= o.fig; 
                setpixelposition(fig_handle(gconi),[1 1 1300 450]);
            end
            con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
            con_report=osl_report_print_figs(con_report);
            
            %%%%%%%%%%%%%%%%%%%
            %% lower level cope plots
            if(oat.group_level.report.show_lower_level_cope_maps)
                fig_handle=[];
                fig_name={};
                fig_title={};
                    
                for subi=1:length(oat.group_level.subjects_to_do),
                    sub=oat.group_level.subjects_to_do(subi);
                    
                    %%%%%%%%%%%%%%%%%%%
                    %% output required niis
                    % output required niis
                    S2=[];
                    S2.oat=oat;
                    S2.time_range=[group_level_results.times(time_ind_max) group_level_results.times(time_ind_max)];
                    S2.stats_fname=oat.subject_level.results_fnames{sub};
                    S2.first_level_contrasts=[con]; % list of first level contrasts to output
                    resamp_gridstep=2;
                    S2.resamp_gridstep=resamp_gridstep;%oat.source_recon.gridstep;
                    S2.freq_bin=freq_ind_max;  
                    S2.stats_dir=report_dir;
                    
                    [statsdir,times,count]=oat_save_nii_stats(S2);

                    fig_name{subi}=['lower_level_cope_at_maxt_smap_c' num2str(con) '_sub' num2str(sub)];  
                    fig_title{subi}=['Lower level cope' num2str(con) ' for sub' num2str(sub) ', at t=' num2str(group_level_results.times(time_ind_max)) 'secs, f=' mat2str(group_level_results.frequency_ranges(freq_ind_max,:)) 'Hz'];  
                    S2=[];                
                    S2.percrange=[96 99.9];
                    %S2.percrange=[70 99.9];
                    S2.mni_coord=group_level_results.mni_coords(vox_ind_max,:);
                    S2.title='';
                    S2.fname=[statsdir '/cope' num2str(con) '_' num2str(resamp_gridstep) 'mm.nii.gz'];                              
                    o = ortho_overlay_act( S2 );
                    fig_handle(subi)= o.fig; 
                    setpixelposition(fig_handle(coni),[1 1 1300 450]);
                end
                con_report=osl_report_set_figs(con_report,fig_name,fig_handle,fig_title);        
                con_report=osl_report_print_figs(con_report);

            end

        end % if ~oat.group_level.space_average,
        
    end % if is_sensor_space 
    
    %%%%%%%%%%%%%%%%%%%
    %% generate source recon web report
    con_report=osl_report_write(con_report, group_level_report);
    group_level_report=osl_report_add_sub_report(group_level_report, con_report);

end % for cc=1:length(first_level_cons_to_do), 

%%%%%%%%%%%%%%%%%%%
%% generate source recon web report
group_level_report=osl_report_write(group_level_report, oat.results.report);        
disp(['To view group level report, point your browser to <a href="' group_level_report.html_fname '">' group_level_report.html_fname '</a>']);

