function topos = sensors_topoplot(D,comp,modality,do_plot)

% topos = sensors_topoplot(D,comp,modality)
%
% private function to produce sensor topoplot
% can take a cell array of modalities and plots them consecutively
% RB 2017, adapted from MWW 2016

cfg  = [];
data = [];

if nargin<4
    do_plot=0;
end

comp(D.badchannels,:) = 0;

nmodalities=length(modality);

for modidx=1:nmodalities
    
    subplot(1,nmodalities,modidx)
    
    cmodality=modality{modidx};
    
    comp2view = comp(indchantype(D,cmodality),:);
    
    if (strcmp(cmodality,'MEGPLANAR')) % Average gradiometers
        comp2view = sqrt(comp2view(1:2:end,:).^2 + comp2view(2:2:end,:).^2);
    end
    
    if (strcmp(cmodality,'MEGMAG'))
        cfg.channel     = {'MEGMAG'};
        cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');
    elseif (strcmp(cmodality,'MEGPLANAR'))
        cfg.channel     = {'MEGMAG'};
        cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');
    elseif (strcmp(cmodality,'MEGGRAD'))
        cfg.channel     = {'MEG'};
        cfg.layout      = fullfile(osldir,'layouts','CTF275.lay');
    elseif strcmp(cmodality, 'MEG')
        cfg.channel = {'MEG'};
        cfg.layout  = fullfile(osldir,'layouts','4D248.lay');
    elseif (strcmp(cmodality,'EEG')),
        warning('EEG not currently supported, using development EEG layout');
        cfg.channel = {'EEG'};
        cfg.layout  = fullfile(osldir,'layouts','EEG60.lay');
    else
        error('Unsupported modality');
    end
    
    data.dimord    = 'chan_comp';
    data.topo      = comp2view;
    data.topolabel = D.chanlabels(indchantype(D,cfg.channel));
    data.time      = {1};
    
    %cfg = rmfield(cfg,'channel');
    %cfg.component   = 1:size(comp,2);
    cfg.interactive = 'no';
    cfg.comment     = cmodality;
    cfg.commentpos     = 'title';

    cfg.title       = cmodality;
    
    %cfg.layout = ft_prepare_layout(cfg);
    
    if do_plot
        ft_topoplotER(cfg,data);
    else
        tmp_fig = figure('visible','off');
        [~] = evalc('ft_topoplotER(cfg,data);');
        topos = handle2struct(get(gcf,'children'));
        topos = topos(end:-1:1); % handles are LIFO
        close(tmp_fig)
    end
    
end

end



