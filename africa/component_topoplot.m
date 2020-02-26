function topos = component_topoplot(D,comp,modality,do_plot,layout)

% topos = component_topoplot(D,comp,modality)
%
% private function to produce component topoplot used by africa
%
% MWW 2016

if nargin<5
    layout = [];
end

if nargin<4
    do_plot=0;
end


cfg  = [];
data = [];


comp(D.badchannels,:) = 0;
comp2view = comp(indchantype(D,modality),:);

if (strcmp(modality,'MEGPLANAR')) % Average gradiometers
    comp2view = sqrt(comp2view(1:2:end,:).^2 + comp2view(2:2:end,:).^2);
end

if strcmp(modality,'MEGMAG') && strcmp(D.sensors('MEG').type,'neuromag306')
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');

elseif strcmp(modality,'MEGPLANAR') && strcmp(D.sensors('MEG').type,'neuromag306')
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');

elseif strcmp(modality,'MEGGRAD') && strcmp(D.sensors('MEG').type,'ctf275')
    cfg.channel     = {'MEG'};
    cfg.layout      = fullfile(osldir,'layouts','CTF275.lay');

elseif strcmp(modality,'MEGMAG') && strcmp(D.sensors('MEG').type,'bti248')
    cfg.channel = {'MEGMAG'};
    cfg.layout  = fullfile(osldir, 'layouts', '4D248.lay');

elseif (strcmp(modality,'EEG'))
    if isempty(layout)
        warning('EEG layout is being automatically generated from sensor info');
        S = [];
        S.elec = D.sensors('EEG');
        cfg.layout  = ft_prepare_layout( S );
    else
        cfg.layout = layout;
    end
    cfg.channel = {'EEG'};
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



