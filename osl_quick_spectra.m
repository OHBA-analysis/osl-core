function fig = osl_quick_spectra(D,bands,modality)
%% function fig = osl_quick_spectra(D,bands,modality)
%
% Make a quick summary of sensor topologies as a function of frequency
%
% Inputs
% D: meeg object containing sensorspace data
% bands: optional cell array of frequency bands to use
%        defaults to = {[0 3], [3 8], [8 13], [13 30], [30 100], [100 200]}
% modality: optional specification of sensor modality to plot
%        defaults to 'MEGPLANAR'

if nargin < 3 || isempty(modality)
    modality = 'MEGPLANAR';
end

if nargin < 2 || isempty(bands)
    bands = {[0 3], [3 7], [7 13], [13 30], [30 100], [100 200]};
end

chaninds = D.indchantype('MEGPLANAR');
[pxx,f] = pwelch( zscore(D(chaninds,:,1)'),D.fsample,[],[],D.fsample );

font_size = 12;

fig = figure('Position',[100 100 512*2 768]);
subplot(2,2,1);
loglog(f,pxx);grid on;
xlim([1,D.fsample/2]);
yl = ylim;
ylim([mean(mean(pxx(end-50:end,:))) yl(2)])
tks = [1 5 10 20 30 40 50 100 200];
set(gca,'XTick',tks,'XTickLabels',tks,'fontsize',font_size)
xlabel('Frequency (Hz)','fontsize',font_size);
ylabel('Frequency (Hz)','fontsize',font_size);
title('Welch log PSD','fontsize',font_size)
subplot(2,2,2)
cc = corrcoef(pxx');
contourf(log(f(2:end)),log(f(2:end)),cc(2:end,2:end),96,'linestyle','none');
set(gca,'XTick',log(tks),'XTickLabels',tks,'fontsize',font_size)
set(gca,'YTick',log(tks),'YTickLabels',tks,'fontsize',font_size)
xlabel('Frequency (Hz)','fontsize',font_size);
ylabel('Frequency (Hz)','fontsize',font_size);
title('Spatial correlation over frequency','fontsize',font_size)
grid on;
for ii = 1:length(bands)
    subplot(2,length(bands),ii+length(bands));
    freq_inds = f>bands{ii}(1) & f<bands{ii}(2);
    quick_topo(D,squeeze(mean(pxx(freq_inds,:),1))',...
                {modality},[num2str(bands{ii}) ' Hz']);
end

end

function quick_topo(D,comp,modality,sensor_title)


cfg  = [];

if strcmp(modality,'MEGMAG') && strcmp(D.sensors('MEG').type,'neuromag306')
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');

elseif strcmp(modality,'MEGPLANAR') && strcmp(D.sensors('MEG').type,'neuromag306')
    % Combine grads and pretend we have mags - easier than swapping labels
    if (strcmp(modality,'MEGPLANAR')) % Average gradiometers
        comp = sqrt(comp(1:2:end,:).^2 + comp(2:2:end,:).^2);
    end
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = fullfile(osldir,'layouts','neuromag306mag.lay');

elseif strcmp(modality,'MEGGRAD') && strcmp(D.sensors('MEG').type,'ctf275')
    cfg.channel     = {'MEG'};
    cfg.layout      = fullfile(osldir,'layouts','CTF275.lay');

elseif strcmp(modality,'MEGMAG') && strcmp(D.sensors('MEG').type,'bti248')
    cfg.channel = {'MEGMAG'};
    cfg.layout  = fullfile(osldir, 'layouts', '4D248.lay');

elseif (strcmp(modality,'EEG'))
    warning('EEG not currently supported, using development EEG layout');
    cfg.channel = {'EEG'};
    cfg.layout  = fullfile(osldir, 'layouts', 'EEG60.lay');
else
    error('Unsupported modality');
end

data = [];
data.dimord    = 'chan_comp';
data.topo      = comp;
data.topolabel = D.chanlabels(indchantype(D,cfg.channel));
data.time      = {1};

cfg = rmfield(cfg,'channel');
cfg.component   = 1:size(comp,2);
cfg.interactive = 'no';
cfg.comment     = 'no';
cfg.title       = sensor_title;
cfg.colorbar    = 'South';

ft_topoplotIC(cfg,data);
end
