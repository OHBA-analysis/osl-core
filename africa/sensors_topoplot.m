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
    
    if (strcmp(cmodality,'MEGPLANAR')) % Average gradiometers (because there are two of them at each position)
        comp2view = sqrt(comp2view(1:2:end,:).^2 + comp2view(2:2:end,:).^2);
    end
    
    if (strcmp(cmodality,'MEGMAG'))
        %cfg.channel     = {'MEGMAG'};
                cfg.channel     = {'all'};

        cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');
    elseif (strcmp(cmodality,'MEGPLANAR'))
        %cfg.channel     = {'MEGMAG'};
        cfg.channel     = {'all'};

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
    
    % gett
    if (strcmp(cmodality,'MEGMAG')) || (strcmp(cmodality,'MEGPLANAR'))
    data.topolabel = D.chanlabels(indchantype(D,'MEGMAG'));
    else
            data.topolabel = D.chanlabels(indchantype(D,cfg.channel));
    end
    data.time      = {1};

    %if (strcmp(cmodality,'MEGPLANAR')) % Plot every second gradiometer (because there are two of them at each position)
    %    data.topolabel=data.topolabel(1:2:end);
    %end
    
    %cfg = rmfield(cfg,'channel');
    %cfg.component   = 1:size(comp,2);
    cfg.interactive = 'no';
    cfg.comment     = cmodality;
    cfg.commentpos     = 'title';

    cfg.title       = cmodality;
    cfg.skipscale = 'yes'
    cfg.skipcomnt = 'yes'
    
    
    cfg.layout = ft_prepare_layout(cfg);
    
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



