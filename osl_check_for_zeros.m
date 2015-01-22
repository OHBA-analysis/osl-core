function D = osl_check_for_zeros(fname)

D=spm_eeg_load(fname);
[~,chan_inds] = setdiff(find(strcmp(D.chantype, 'MEGMAG') |strcmp(D.chantype, 'MEGPLANAR')| strcmp(D.chantype, 'MEGGRAD')),D.badchannels);
BadEpochs = load_meeg(D);
tmp=sum(D(chan_inds,:,:),1)==0;
if sum(tmp)~=0;
    i=1;eyes_open=1;
    while i<D.nsamples+1;
        if tmp(i)==1 && eyes_open
            if i==1
                BadEpochs{end+1}(1)=D.time(i);
            else
                BadEpochs{end+1}(1)=D.time(i-1);
            end
            eyes_open=0;
        end
        if tmp(i)~=1 && ~eyes_open
            BadEpochs{end}(2)=D.time(i);
            eyes_open=1;
        end
        i=i+1;
    end
    BadEpochs=BadEpochs(cellfun(@numel,BadEpochs)==2);
    D = save_meeg(D,BadEpochs);
end
end
%% Loads Bad epoch information from SPM object
function BadEpochs = load_meeg(D)
BadEpochs = {};
Events = D.events;
for ev = 1:numel(Events)
    if strncmp({Events.type},'artefact',8)
        BadEpochs{end+1}(1) = Events(ev).time;
        BadEpochs{end}(2) = Events(ev).time + Events(ev).duration;
    end
end
end
%% Saves Bad epoch information to SPM object
function D = save_meeg(D,BadEpochs)
Events = struct([]);
% Save bad epochs using method meeg/events
for ev = 1:numel(BadEpochs)
    if numel(BadEpochs{ev} == 2)
        Events(ev).type     = 'artefact_zeros';
        Events(ev).value   = ev;
        Events(ev).time     =  BadEpochs{ev}(1);
        Events(ev).duration = diff(BadEpochs{ev});
    end
end
D = events(D,1,Events);
save(D);
end