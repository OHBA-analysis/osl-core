function [fig_handle, fig_name, fig_title] = oat_plot_vox_stats(S)

% [fig_handle fig_name fig_title] = oat_plot_vox_stats(S)
%
% Plot the copes and t-stats at a specified vox coordinate
%
% S.stats; % oat.first_level or oat.group_level results struct or fname
% (for which you will also need to specify S.oat)
%
% S.vox_coord; % mni coordinate to use if in source space, 
% OR
% S.chanlabel; % if in sensor space
% OR
% neither if stats.cope only contains one "voxel" or "channel"
%
% S.freq_inds; % only plot these freqs
% S.first_level_cons_to_do; % list of first level cons to plot
% S.group_level_cons_to_do; % list of group level cons to plot (if group
% level stats)
%
% MWW 2013

try stats=S.stats; catch, error('S.stats not specified'); end;

if isstr(stats)
    stats=oat_load_results(S.oat,stats);
else
    stats=stats;
end;

try first_level_cons_to_do=S.first_level_cons_to_do; catch, first_level_cons_to_do=1:size(stats.cope,3); end;
try group_level_cons_to_do=S.group_level_cons_to_do; catch, group_level_cons_to_do=1:size(stats.cope,5); end;

try freq_inds=S.freq_inds; catch, freq_inds=1:length(stats.frequencies); end; % only plot these freqs
try vox_coord=S.vox_coord; catch, vox_coord=[]; end;

if stats.level==2,    
    is_sensor_space=strcmp(stats.recon_method,'none');
else
    is_sensor_space=strcmp(stats.source_recon.method,'none');
end;

if(stats.level==2),
    first_level_cons_to_do=intersect(first_level_cons_to_do,stats.current_level.first_level_contrasts_to_do);
    
    %if(length(group_level_cons_to_do)~=length(stats.current_level_contrast_name)),
    %    error('Number of group_level_cons_to_do not compatible with number of contrasts');
    %end;
    
    first_level_contrast_name=stats.lower_level_contrast_name;
    group_level_contrast_name=stats.current_level_contrast_name;
else
    first_level_contrast_name=stats.first_level_contrast_name;
    group_level_contrast_name={};    
end;

times=stats.times;
frequencies=stats.frequencies(freq_inds);

if size(stats.cope,1)==1 && isempty(vox_coord),    
    vox_ind=1;   
    spatial_loc_str='given loc';
else
    if ~is_sensor_space,        
            
        if isempty(vox_coord),            
            error('S.vox_coord needs to be specified');
        end;
        vox_coord=S.vox_coord;

        % Find the nearest index of the beamformed voxels to the specified mni_coord 
        [vox_ind,vec,dist]=nearest_vec(stats.mni_coords,vox_coord);
        spatial_loc_str=mat2str(stats.mni_coords(vox_ind,:),4);
        
    else
        if isempty(S.chanlabel),
            error('S.chanlabel needs to be specified');
        end;
        chanlabel=S.chanlabel;

        % Find the ind for the specified chan label       
        vox_ind = find(strcmp(stats.chanlabels, chanlabel));
        
        spatial_loc_str=chanlabel;
       
    end;
        
end;

ntpts = length(times);
nfreqs = length(frequencies);

if(ntpts==1)
    error('Can not do time series plots with only one timepoint');    
end;

cols={'r','g','b','y','m','k'};
fig_handle=[];
fig_name={};
fig_title={};

% use first contrast to choose the voxel and tpts of interest
con=first_level_cons_to_do(1);
% use first group contrast as one to choose coords
gcon=group_level_cons_to_do(1);

%%%%%%%
%% do plots of first level cons shown together for each group contrast
for gconi=1:length(group_level_cons_to_do),
    
    gcon=group_level_cons_to_do(gconi);
    
    fig_handle(gconi)=sfigure;
    
    if(nfreqs==1),
        set(fig_handle(gconi),'Position',[1 1 1300 450]);   
    else
        set(fig_handle(gconi),'Position',get(fig_handle(gconi),'Position')*2);
    end;
    
    if(stats.level==2),
        fig_name{gconi}=['stats_tc_gc' num2str(gcon) '_' group_level_contrast_name{gcon} ' at ' spatial_loc_str];
        fig_title{gconi}=['Stats for gc' num2str(gcon) ', ' group_level_contrast_name{gcon} ' at ' spatial_loc_str];
        
    else
        fig_name{gconi}='stats_tc';
        fig_title{gconi}=['Stats at ' spatial_loc_str];
    end;
    
    leg={};

    max_cons=6;
    if(length(first_level_cons_to_do)>max_cons)
        disp(['Only plotting the first ' num2str(max_cons) ' contrasts']);
        first_level_cons_to_do=first_level_cons_to_do(1:max_cons);
    end;
    
    if(nfreqs>1),        
        
        %%%%%%%%%%%%%
        % plot TF images              
        for coni=1:length(first_level_cons_to_do),
            con=first_level_cons_to_do(coni);    

            cope=permute(stats.cope(vox_ind,:,con,freq_inds,gcon),[4 2 1 3 5]);
            stdcope=permute(stats.stdcope(vox_ind,:,con,freq_inds,gcon),[4 2 1 3 5]);
                
            % plot copes        
            subplot(length(first_level_cons_to_do),2,(coni-1)*2+1);                        
            imagesc(times,frequencies,cope);
            axis xy; colorbar;
            title([spatial_loc_str ', cope' num2str(con) ' [' first_level_contrast_name{con} ']']);
            plot4paper('time (s)','freq (Hz)');
            
            % plot 1-tailed tstats
            tstat=cope./stdcope;
            subplot(length(first_level_cons_to_do),2,(coni-1)*2+2);                        
            imagesc(times,frequencies,tstat);
            axis xy; colorbar;
            title([spatial_loc_str ', tstat' num2str(con) ' [' first_level_contrast_name{con} ']']);
            plot4paper('time (s)','freq (Hz)');
        end;       

    else,
       
        % do time series plots
        for coni=1:length(first_level_cons_to_do),
            con=first_level_cons_to_do(coni);

            cope=permute(stats.cope(vox_ind,:,con,freq_inds,gcon),[2 1 3 4 5]);
            stdcope=permute(stats.stdcope(vox_ind,:,con,freq_inds,gcon),[2 1 3 4 5]);
        
            % plot copes
            subplot(1,2,1);
            hold on;
            errorbar(times,cope,stdcope,cols{coni});

            % plot 1-tailed tstats
            tstat=cope./stdcope;
            subplot(1,2,2);
            hold on;
            plot(times, tstat, cols{coni},'LineWidth',2); hold on;   

            leg{length(leg)+1}=first_level_contrast_name{con};

        end;

        subplot(1,2,1);
        a=axis; axis([times(1) times(end) a(3) a(4)]);    
        plot4paper('time (secs)','cope'); 
        legend(leg,'Location','NorthWest');         
        titstr=[spatial_loc_str];
        if(stats.level==2)
           titstr=[titstr, ', gc' num2str(gcon) ', ' group_level_contrast_name{gcon}];
        end;
        title(titstr);
        
        subplot(1,2,2);
        a=axis; axis([times(1) times(end) a(3) a(4)]);    
        plot4paper('time (secs)','t-stat'); 
        legend(leg,'Location','NorthWest'); 
        titstr=[spatial_loc_str];
        if(stats.level==2)
           titstr=[titstr, ', gc' num2str(gcon) ', ' group_level_contrast_name{gcon}];
        end;
        title(titstr);
        
    end;
    
end;

figindstart=length(fig_handle);

if(stats.level==2 && length(group_level_cons_to_do)>1)

    %%%%%%%
    %% do plots of group level cons shown together for each first level contrast
    for coni=1:length(first_level_cons_to_do),
        con=first_level_cons_to_do(coni);    

        fig_handle(figindstart+coni)=sfigure;

        if(nfreqs==1),
            set(fig_handle(figindstart+coni),'Position',[1 1 1300 450]);   
        else
            set(fig_handle(figindstart+coni),'Position',get(fig_handle(coni),'Position')*2);
        end;

        fig_name{figindstart+coni}=['stats_tc_c' num2str(con)];
        fig_title{figindstart+coni}=['Stats for c' num2str(con) ' [' first_level_contrast_name{con} ']'];

        leg={};

        max_cons=6;
        if(length(first_level_cons_to_do)>max_cons)
            disp(['Only plotting the first ' num2str(max_cons) ' contrasts']);
            first_level_cons_to_do=first_level_cons_to_do(1:max_cons);
        end;

        if(nfreqs>1),        

            %%%%%%%%%%%%%
            % plot TF images              
            for gconi=1:length(group_level_cons_to_do),    
                gcon=group_level_cons_to_do(gconi);

                cope=permute(stats.cope(vox_ind,:,con,freq_inds,gcon),[4 2 1 3 5]);
                stdcope=permute(stats.stdcope(vox_ind,:,con,freq_inds,gcon),[4 2 1 3 5]);

                % plot copes        
                subplot(length(first_level_cons_to_do),2,(coni-1)*2+1);                        
                imagesc(times,frequencies,cope);
                axis xy; colorbar;
                title(['cope gc' num2str(gcon) ', ' group_level_contrast_name{gcon}]);
                plot4paper('time (s)','freq (Hz)');

                % plot 1-tailed tstats
                tstat=cope./stdcope;
                subplot(length(first_level_cons_to_do),2,(coni-1)*2+2);                        
                imagesc(times,frequencies,tstat);
                axis xy; colorbar;
                title(['tstat gc' num2str(gcon) ', ' group_level_contrast_name{gcon}]);
                plot4paper('time (s)','freq (Hz)');
            end;       

        else,

            % do time series plots
            for gconi=1:length(group_level_cons_to_do),    
                gcon=group_level_cons_to_do(gconi);

                cope=permute(stats.cope(vox_ind,:,con,freq_inds,gcon),[2 1 3 4 5]);
                stdcope=permute(stats.stdcope(vox_ind,:,con,freq_inds,gcon),[2 1 3 4 5]);

                % plot copes
                subplot(1,2,1);
                hold on;
                errorbar(times,cope,stdcope,cols{gconi});

                % plot 1-tailed tstats
                tstat=cope./stdcope;
                subplot(1,2,2);
                hold on;
                plot(times, tstat, cols{gconi},'LineWidth',2); hold on;   

                leg{length(leg)+1}=group_level_contrast_name{gcon};

            end;

            subplot(1,2,1);
            a=axis; axis([times(1) times(end) a(3) a(4)]);    
            plot4paper('time (secs)','cope'); 
            legend(leg,'Location','NorthWest'); 
            
            titstr=[spatial_loc_str];
            if(stats.level==2)
               titstr=[titstr, ', c' num2str(con) ' [' first_level_contrast_name{con} ']'];
            end;
            title(titstr);

            subplot(1,2,2);
            a=axis; axis([times(1) times(end) a(3) a(4)]);    
            plot4paper('time (secs)','t-stat'); 
            legend(leg,'Location','NorthWest'); 
            titstr=[spatial_loc_str];
            if(stats.level==2)
               titstr=[titstr, ', c' num2str(con) ' [' first_level_contrast_name{con} ']'];
            end;
            title(titstr);

        end;

    end;

end
