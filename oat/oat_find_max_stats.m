function [vox_ind_max time_ind_max freq_ind_max stats max_stat] = oat_find_max_stats( S )

% [vox_ind_max time_ind_max freq_ind_max] = oat_find_max_stats(S)
%
% Find max stat location in space, time and freq
%
% S.stats; % oat.first_level or oat.group_level results struct or fname
% (for which you will also need to specify S.oat)
% S.first_level_cons_to_do; % uses first one in list to find stats
% S.group_level_cons_to_do; % uses first one in list to find stats
% S.max_method % method to use, e.g. 'max_abs_tstat'
% S.modality % Only set if in sensor space, e.g. MEGMAG or MEGPLANAR
% S.time_range % used to restrict max method search
% S.freq_range % used to restrict max method search
%
% MWW 2013

try, stats=S.stats; catch error('S.stats not specified'); end;

if isstr(stats)
    stats=oat_load_results(S.oat,stats);
else
    stats=stats;
end;

try, max_method=S.max_method; catch max_method='max_abs_tstat'; end;

try, first_level_con=S.first_level_cons_to_do(1); catch first_level_con=1; end;
try, group_level_con=S.group_level_cons_to_do(1); catch group_level_con=1; end;

try, time_range=S.time_range; catch time_range=[]; end; % used to restrict max method search
try, freq_range=S.freq_range; catch freq_range=[]; end;  % used to restrict max method search

times=stats.times;
frequencies=stats.frequencies;
ntpts = length(times);
nfreqs = length(frequencies); 

if isempty(time_range),
    time_range=[0 stats.times(end)];    
end;

if ntpts==1
    time_range=[stats.times(1) stats.times(end)];
end;

if isempty(freq_range),
    freq_range=[stats.frequencies(1) stats.frequencies(end)];
end;

if stats.level==2,    
    is_sensor_space=strcmp(stats.recon_method,'none');
else
    is_sensor_space=strcmp(stats.source_recon.method,'none');
end;
      
if is_sensor_space,
    if size(S.stats.cope,1)==1,
        % data been space averaged
        vox_inds=1;
    else
        
        try, modality=S.modality; catch error('Must set S.modality'); end;

        if strcmp(modality,'MEGPLANAR'), modality='MEGCMB'; end;
        if strcmp(modality,'MEG'), modality='MEGGRAD'; end;
        
        vox_inds=find(strcmp(stats.chantype,modality));
    end;
else
    vox_inds=1:size(stats.cope,1);
end;

con=first_level_con;
gcon=group_level_con;

time_inds=find(stats.times>=time_range(1) & stats.times<=time_range(2));
if(stats.frequencies(1)==-1),
    freq_inds=1;
else
    freq_inds=find(stats.frequencies>=freq_range(1) & stats.frequencies<=freq_range(2));
end;

cope=permute(stats.cope(vox_inds,time_inds,con,freq_inds,gcon),[1 2 4 3 5]);
stdcope=permute(stats.stdcope(vox_inds,time_inds,con,freq_inds,gcon),[1 2 4 3 5]);
tstats=cope./stdcope; %nvox x ntime x nfreq

switch max_method,

    case 'max_abs_tstat',            
        tstats=abs(tstats);
        %[vox_ind_max time_ind_max]=max2d(tstats);            
        [vox_ind_max time_ind_max freq_ind_max max_stat]=max3d(tstats);                        
        time_ind_max=time_inds(time_ind_max);
        freq_ind_max=freq_inds(freq_ind_max);
        vox_ind_max=vox_inds(vox_ind_max);
        
    case 'max_tstat',
        %[vox_ind_max time_ind_max]=max2d(tstats); 
        [vox_ind_max time_ind_max freq_ind_max max_stat]=max3d(tstats);                        
        time_ind_max=time_inds(time_ind_max);            
        freq_ind_max=freq_inds(freq_ind_max); 
        vox_ind_max=vox_inds(vox_ind_max);
        
    otherwise
        error('Unknown max_method method');

end;

end

