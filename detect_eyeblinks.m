function D = detect_eyeblinks(S)

% D = detect_eyeblinks(S)
%
% detects eyeblinks in spm continuous data file
%
% inputs: S.D      = spm eeg object
%          .method = 'Events' or 'Topo' (topography) or 'Both' (generally
%                     use 'Events': 'Topo' is very coarse,
%                     better to use PCA)
%          .doplots= plot detected eyeblinks and the average
%                    waveform
%          .sd_thr = threshold to reject things that look like
%                    eye-blinks but probably aren't...

% LH 150709

if nargin == 0
    S = [];
end

if ischar(S)
    D = S;
    S = [];
    S.D = D;
        method = spm_input('What to return','+1', 'Events|Topography|Both', strvcat('events', 'topo', 'both'));
else
    try 
        D = S.D;
    catch
        D = spm_select(1, '\.mat$', 'Select EEG mat file');
        S.D = D;
    end
    try 
        method = S.method; 
    catch 
        method = spm_input('What to return','+1', 'Events|Topography|Both', strvcat('events', 'topo', 'both'));
    end
    try 
        doplots = S.doplots; 
    catch
        doplots = 0;
    end
    try
        sd_thr = S.sd_thr;
    catch
        sd_thr = 3;
    end
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

%% get EOG data

eog = find(strcmp(D.chantype,'EOG')); %identify EOG channel
if isempty(eog)
    error('Couldn''t find eog channel');
end
if length(eog)~=1
    error('More than one EOG channel - not currently supported')
end

eog_data = D(eog,:,:);

%% filter data at 1-15Hz (eyeblink duration typically 100-300ms) and demean
lp = 2*1./D.fsample;
hp = 2*15./D.fsample;
h1=fir1(1001,[lp hp]);
%to view filter properties can use: fvtool(h1,1,'Fs',D.fsample)
eog_filt = demean(filtfilt(h1,1,eog_data));

%% find eye-movements

%sd_eeg = std(eog_filt); %MW

%MW:
sd_eeg=(percentile(eog_filt,50)-percentile(eog_filt,15)+percentile(eog_filt,85)-percentile(eog_filt,50))/2; %MW

em_lthresh = sd_thr*sd_eeg;
em_uthresh = sd_thr*10*sd_eeg;
% END MW

%% find 'spikes' (putative eyeblinks):

eblength = round(D.fsample/5); %length of eyeblink(200 ms) in samples;
spikes = [];
for i = eblength:length(eog_filt)-eblength;
    if abs(eog_filt(i))>em_lthresh && ... %bigger than threshold
        abs(eog_filt(i))<em_uthresh && ...    
       all(abs(eog_filt(i))>=abs(eog_filt(i-eblength+1:i+eblength))); %biggest in surrounding 400ms
        spikes = [spikes i];
    end
end

if isempty(spikes)
    error('No eye-blinks detected by algorithm. Try a lower threshold.')
end

% MW:
num_eb_per_min=60*length(spikes)/(D.time(end)-D.time(1));
disp([num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm'])
if (num_eb_per_min<5)
    error(['Only ' num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm. Try a lower threshold.'])
end
if (num_eb_per_min>60)
    error(['As many as ' num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm. Try a higher threshold.'])
end
% END MW

spikemat = zeros(eblength*2,length(spikes));
for i = 1:length(spikes)
    spikemat(:,i) = eog_filt(spikes(i)-eblength+1:spikes(i)+eblength);
end

%reject spikes whose peak is not within 0.5 s.d. of the mean (gets rid of most artefacts
%    etc. not removed by filtering):
mn_spike = mean(spikemat(eblength,:));
sd_spike = 0.5*std(spikemat(eblength,:));
spikes(spikemat(eblength,:)>mn_spike+sd_spike | ...
       spikemat(eblength,:)<mn_spike-sd_spike) = [];
spikemat(:,find(spikemat(eblength,:)>mn_spike+sd_spike | ...
       spikemat(eblength,:)<mn_spike-sd_spike)) = [];

disp(['Number of putative eyeblinks detected: ' num2str(length(spikes))]);
          
%% plot if requested
if doplots == 1
    sfigure;
    plot(spikes,ones(length(spikes),1)*5*sd_eeg,'r.'); 
    hold on;
    plot(eog_filt);
    print('-dpng',[S.fname '_eb_spikes']);
    
    sfigure; hold on;
    plot(abs(spikemat));plot(mean(abs(spikemat),2),'Color','k','LineWidth',4);
    print('-dpng',[S.fname '_eb_mean_timecourse']);
    
end

%% append spikes to list of events
evs = D.events;
for i = 1:length(spikes);
    evs(end+1).type = 'eyeblink';
    evs(end).time = spikes(i)/D.fsample + D.timeonset;
    evs(end).duration = [];
    evs(end).value = 1;
end
if strcmp(lower(method),'events')||strcmp(lower(method),'both')
    D = events(D,1,evs);
end

%% compute topography at peak of eyeblink and add as spatial confound
if strcmp(lower(method),'topo')||strcmp(lower(method),'both')
    chanlist = D.meegchannels;
    spikemeegdat = D(chanlist,spikes,1);
    spiketop = mean(spikemeegdat,2);
    newsconfound = [];

    newsconfound.label = D.chanlabels(chanlist);
    newsconfound.coeff = spiketop;
    newsconfound.bad = zeros(length(chanlist),1);

    D = sconfounds(D,newsconfound);
end

%% plot if requested
if strcmp(lower(method),'topo')||strcmp(lower(method),'both')
    if doplots == 1
        view_spatial_confounds(D);
    end
end

    





    