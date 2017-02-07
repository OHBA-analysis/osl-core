function [cfg, data, fig_handle]=oat_stats_multiplotTFR( S )

% [cfg, dats, fig_handle]=oat_stats_multiplotTFR(S)
%
% calls ft_multiplotTFR passing on S.cfg settings to display the tstats for contrast S.contrast for
% OAT S.stats
%
% Input:
% S.oat
% S.stats_fname 
% S.modality, e.g. 'MEGPLANAR'
% S.first_level_contrast. Specify single first level contrast index
% S.group_level_contrast. Specify single group contrast index (if group level stats)
% S.cfg 
% S.view_cope. Binary flag, if 0 (default) then tstats are viewed
%
% MWW 2012

OSLDIR = getenv('OSLDIR');

try, oat=S.oat; catch, error('No S.oat supplied'); end;
try, stats_fname=S.stats_fname; catch, error('No S.stats_fname supplied'); end;
try, modality=S.modality; catch, modality='MEGPLANAR'; end;
try, contrast=S.first_level_contrast; catch, error('No S.first_level_contrast supplied'); end; % first level contrast - there can be multiple numbers of these
try, group_contrast=S.group_level_contrast; catch, group_contrast=1; end; % specify group contrast if group level stats
try, view_cope=S.view_cope; catch, view_cope=0; end; 
try, cfg=S.cfg; catch, cfg=[]; end;
try, do_plots=S.do_plots; catch, do_plots=1; end;
try chans=S.chans; catch, chans = [];end

try, D_tstat=S.D; catch
    
    if isstr(stats_fname)
        stats=oat_load_results(oat,stats_fname);
    else
        stats=stats_fname;
    end;
    
    S4=[];
    S4.oat=oat;
    S4.stats=stats;
    S4.write_copes=1;
    
    [D_tstat, D_cope]=oat_save_spm_stats(S4);    

    if(stats.level==2),
        D_tstat=D_tstat{group_contrast};
        D_cope=D_cope{group_contrast};
    end;
end;

if(length(contrast)>1)
    error('Can only specify one first level contrast index');
end;

cope=stats.cope;
stdcope=stats.stdcope;

if(stats.level==2),
    cope = cope(:,:,:,:,group_contrast);
    stdcope = stdcope(:,:,:,:,group_contrast);
end;

if isfield(S,'time_range');
    time_inds = find(stats.times > S.time_range(1) & stats.times < S.time_range(2));
else
    time_inds = 1:size(cope,2);
end

data = [];
data.time =  stats.times(time_inds);
data.freq= stats.frequencies;
data.dimord='chan_freq_time';

if(strcmp(modality,'MEGMAG')),
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
    chanind = strmatch('MEGMAG', D_tstat.chantype);
    
    data.label=D_tstat.chanlabels(chanind);
    for ind=1:length(chanind),
        if(~view_cope)
            data.powspctrm(ind,:,:)=permute(cope(chanind(ind),time_inds,contrast,:)./stdcope(chanind(ind),time_inds,contrast,:),[1 4 2 3]);
        else
            data.powspctrm(ind,:,:)=permute(cope(chanind(ind),time_inds,contrast,:),[1 4 2 3]);
        end;
    end;
    
elseif (strcmp(modality,'MEGPLANAR')),
    
    cfg.layout      = [OSLDIR '/layouts/neuromag306cmb.lay'];
    chanind = strmatch('MEGCMB', D_tstat.chantype);
    
    data.label=D_tstat.chanlabels(chanind);
    for ind=1:length(chanind),
        if(~view_cope)
            data.powspctrm(ind,:,:)=permute(cope(chanind(ind),time_inds,contrast,:)./stdcope(chanind(ind),time_inds,contrast,:),[1 4 2 3]);
        else
            data.powspctrm(ind,:,:)=permute(cope(chanind(ind),time_inds,contrast,:),[1 4 2 3]);
        end;
    end;
    
elseif (strcmp(modality,'MEG')) || (strcmp(modality,'MEGGRAD')),
    
    cfg.layout      = [OSLDIR '/layouts/CTF275.lay'];
    chanind = strmatch('MEGGRAD', D_tstat.chantype);
    
    data.label=D_tstat.chanlabels(chanind);
    
    for ind=1:length(chanind),
        indp=ind;
        if(~view_cope)
            data.powspctrm(ind,:,:)=permute(cope(indp,time_inds,contrast,:)./stdcope(indp,time_inds,contrast,:),[1 4 2 3]);
        else
            data.powspctrm(ind,:,:)=permute(cope(indp,time_inds,contrast,:),[1 4 2 3]);
        end;
    end;
    
elseif(strcmp(modality,'EEG')),  % added by DM
    
    B=sensors(D_tstat,'EEG');
    mat=spm2fieldtrip(D_tstat);
    cfgx=[];
    cfgx.elec.pnt=B.elecpos;
    cfgx.elec.label=B.label;
    lay = ft_prepare_layout(cfgx,mat);
    
    cfg.layout =lay;
    chanind = strmatch('EEG', D_tstat.chantype);
    
    data.label=D_tstat.chanlabels(chanind);
    
    for ind=1:length(chanind),
        indp=ind;
        if(view_cope)
            data.powspctrm(ind,:,:)=permute(cope(indp,time_inds,contrast,:),[1 4 2 3]);
        else
            data.powspctrm(ind,:,:)=permute(cope(indp,time_inds,contrast,:)./stdcope(indp,time_inds,contrast,:),[1 4 2 3]);
        end;
    end;
else
    error('Unknown modality');
end;

data.powspctrm=double(data.powspctrm);

if(~isfield(cfg,'interactive')),
    cfg.interactive = 'yes';
end;
if(~isfield(cfg,'zlim')),
    %cfg.zlim = [-4 4];
end;
if(~isfield(cfg,'xlim')),
    cfg.xlim = [stats.times(1) stats.times(end)];
end;

cfg.comment  = '';

cfg.colorbar = 'yes';


fig_handle=[];
if do_plots
    if ~isempty(chans)
        cfg.channel = chans;
        fig_handle = sfigure; ft_singleplotTFR(cfg,data);
    else
        fig_handle = sfigure; ft_multiplotTFR(cfg,data);
    end
end;


