function [cfg, dats, fig_handle]=oat_stats_multiplotER(S)

% [cfg, dats, fig_handle]=oat_stats_multiplotER(S)
%
% calls ft_multiplotER (passing on cfg settings) to display the tstats for contrast S.contrast for
% OAT S.stats
%
% Input:
% S.oat
% S.stats_fname % can be file name of first_level_results struct, or can directly be the first_level_results struct, from
% running a oat.first_level
% S.modality, e.g. 'MEGPLANAR'
% S.first_level_contrast. First level contrast indices - there can be multiple numbers of these
% S.group_level_contrast. Specify single group contrast index (if group level stats)
% S.cfg
% S.view_cope. Binary flag, if 0 (default) then tstats are viewed
%
% MWW 2012

OSLDIR = getenv('OSLDIR');

try, oat=S.oat; catch, error('No S.oat supplied'); end
try, stats_fname=S.stats_fname; catch, error('No S.stats_fname supplied'); end
try, modality=S.modality; catch, modality='MEGPLANAR'; end
try, contrast=S.first_level_contrast; catch, error('No S.first_level_contrast supplied'); end
try, group_contrast=S.group_level_contrast; catch, group_contrast=1; end % specify group contrast if group level stats
try, view_cope=S.view_cope; catch, view_cope=0; end
try, cfg=S.cfg; catch, cfg=[]; end
try, do_plots=S.do_plots; catch, do_plots=1; end
try plot_avg=S.plot_avg; catch, plot_avg=0; end

try, D_tstat=S.D; catch

    if isstr(stats_fname)
        stats=oat_load_results(oat,stats_fname);
    else
        stats=stats_fname;
    end

    S4=[];
    S4.oat=oat;

    S4.stats=stats;
    S4.write_copes=1;

    [D_tstat, D_cope]=oat_save_spm_stats(S4);

    if(stats.level==2),
        D_tstat=D_tstat{group_contrast};
        D_cope=D_cope{group_contrast};
    end
end

cope=stats.cope;
if(stats.level==2),
    cope = cope(:,:,:,:,group_contrast);
end

if(isfield(S,'data')),
    cope=S.data;
    view_cope=1;
    contrast=1;
end

if isfield(S,'time_range');
    time_inds = find(stats.times > S.time_range(1) & stats.times < S.time_range(2));
else
    time_inds = 1:size(cope,2);
end

for cc=1:length(contrast),

    data.time = stats.times(time_inds);
    data.dimord='chan_time';

    if(strcmp(modality,'MEGMAG')),
        cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];

%         cfg2 = [];
%         cfg2.layout = cfg.layout;
%         laygrad = ft_prepare_layout(cfg2);

        chanind         = strmatch('MEGMAG', D_tstat.chantype);

        data.label=D_tstat.chanlabels(chanind);

        ts=zeros(numel(data.label),length(time_inds));
        for ind=1:length(chanind),
            if(view_cope)
                ts(ind,:)=cope(chanind(ind),squeeze(time_inds),contrast(cc));
            else
                ts(ind,:)=cope(chanind(ind),time_inds,contrast(cc))./stats.stdcope(chanind(ind),time_inds,contrast(cc));
            end
        end

    elseif (strcmp(modality,'MEGPLANAR')),

        cfg.layout      = [OSLDIR '/layouts/neuromag306cmb.lay'];
        chanind = strmatch('MEGCMB', D_tstat.chantype);

        data.label=D_tstat.chanlabels(chanind);

        ts=zeros(numel(chanind),length(time_inds));
        for ind=1:length(chanind),
            if(view_cope)
                ts(ind,:)=cope(chanind(ind),time_inds,contrast(cc));
            else
                ts(ind,:)=cope(chanind(ind),time_inds,contrast(cc)) ./ (stats.stdcope(chanind(ind),time_inds,contrast(cc)));
            end
        end

    elseif (strcmp(modality,'MEG')) || (strcmp(modality,'MEGGRAD')),

        cfg.layout      = [OSLDIR '/layouts/CTF275.lay'];
        chanind = strmatch('MEG', D_tstat.chantype);

        data.label=D_tstat.chanlabels(chanind);

        ts=zeros(length(chanind),length(time_inds));
        for ind=1:length(chanind),
            indp=ind;
            if(view_cope)
                ts(ind,:)=cope(indp,time_inds,contrast(cc));
            else
                ts(ind,:)=cope(indp,time_inds,contrast(cc))./stats.stdcope(indp,time_inds,contrast(cc));
            end
        end

    elseif(strcmp(modality,'EEG')),   % added by DM

        B=sensors(D_tstat,'EEG');
        mat=spm2fieldtrip(D_tstat);
        cfgx=[];
        cfgx.elec.pnt=B.elecpos;
        cfgx.elec.label=B.label;
        lay = ft_prepare_layout(cfgx,mat);

        cfg.layout =lay;
        chanind = strmatch('EEG', D_tstat.chantype);

        data.label=D_tstat.chanlabels(chanind);

        ts=zeros(length(chanind),length(time_inds));
        for ind=1:length(chanind),
            if(view_cope)
                ts(ind,:)=cope(ind,time_inds,contrast(cc));
            else
                ts(ind,:)=cope(ind,time_inds,contrast(cc))./stats.stdcope(ind,time_inds,contrast(cc));
            end
        end

    else
        error('Unsupported modality');
    end

    data.avg = ts;
    %data.datatype='timelock';

    dats{cc}=data;
end

if(~isfield(cfg,'interactive')),
    cfg.interactive = 'yes';
end
if(~isfield(cfg,'xlim')),
    cfg.xlim        = [stats.times(time_inds(1)) stats.times(time_inds(end))];
end
if(~isfield(cfg,'comment')),
    cfg.comment     = '';
end
cfg.parameter='avg';

fig_handle=[];

if do_plots,
    if plot_avg == 1  || length(dats{1}.time) == 1
        fig_handle=sfigure; ft_topoplotER(cfg,dats{:});
    else
        fig_handle=sfigure; ft_multiplotER(cfg,dats{:});
    end
end
