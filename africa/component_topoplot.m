function topos = component_topoplot(D,comp,modality,do_plot)

% topos = component_topoplot(D,comp,modality)
%
% private function to produce component topoplot used by africa
%
% MWW 2016

cfg  = [];
data = [];
 
if nargin<4
    do_plot=0;
end

global OSLDIR

comp(D.badchannels,:) = 0;
comp2view = comp(indchantype(D,modality),:);

if (strcmp(modality,'MEGPLANAR')) % Average gradiometers
    comp2view = sqrt(comp2view(1:2:end,:).^2 + comp2view(2:2:end,:).^2);
end

if (strcmp(modality,'MEGMAG'))
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
elseif (strcmp(modality,'MEGPLANAR'))
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
elseif (strcmp(modality,'MEGGRAD'))
    cfg.channel     = {'MEG'};
    cfg.layout      = [OSLDIR '/layouts/CTF275.lay'];
elseif strcmp(modality, 'MEG')
    cfg.channel = {'MEG'};
    cfg.layout  = fullfile(OSLDIR, 'layouts', '4D248.lay');
elseif (strcmp(modality,'EEG')),
    warning('EEG not currently supported, using development EEG layout');
    cfg.channel = {'EEG'};
    cfg.layout  = [OSLDIR '/layouts/EEG60.lay'];
else
    error('Unsupported modality');
end

data.dimord    = 'chan_comp';
data.topo      = comp2view;
data.topolabel = D.chanlabels(indchantype(D,cfg.channel));
data.time      = {1};

cfg = rmfield(cfg,'channel');
cfg.component   = 1:size(comp,2);
cfg.interactive = 'no';
cfg.comment     = 'no';
cfg.title       = modality;

%cfg.layout = ft_prepare_layout(cfg);

if do_plot
    ft_topoplotIC(cfg,data);
else    
    tmp_fig = figure('visible','off');
    [~] = evalc('ft_topoplotIC(cfg,data);');
    topos = handle2struct(get(gcf,'children'));
    topos = topos(end:-1:1); % handles are LIFO
    close(tmp_fig)
end

end



