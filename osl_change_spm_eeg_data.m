function [ D2 ] = osl_change_spm_eeg_data( Sc )

% [ D2 ] = osl_change_spm_eeg_data( Sc )
%
% Clones a passed in Sc.D and fills it with the passed in
% data. Needs:
% Sc.D 
% Sc.newdata (3D matrix of data (channels, samples, trials))
% Sc.newname (filename - should end with suffix .dat)
%
% Optional:
% Sc.time (vector of time points)
% Sc.cond_list
% Sc.frequencies (vector of freqs)
% Sc.modalities ('MEG' or 'EEG')
% Sc.remove_montages (default is yes)
% Sc.chantype (default is 'LFP')
% 
% MWW

try, tmp=Sc.cond_list; catch, Sc.cond_list=Sc.D.conditions; end;
try, modalities=Sc.modalities; catch, modalities=[];  end;
try, Sc.time=Sc.time; catch, Sc.time=Sc.D.time; end;
try, Sc.remove_montages=Sc.remove_montages; catch, Sc.remove_montages=1; end;
try, chantype=Sc.chantype; catch, chantype='LFP';  end;

% remove montages
if Sc.remove_montages
    D=montage(Sc.D,'remove',1:montage(Sc.D,'getnumber'));
else
    D=Sc.D;
end

if(size(Sc.newdata,4)==1)
    %TF    
    D2=clone(D,Sc.newname,[size(Sc.newdata,1),size(Sc.newdata,2),size(Sc.newdata,3)],0);
else
    D2=clone(D,Sc.newname,[size(Sc.newdata,1),size(Sc.newdata,4),size(Sc.newdata,2),size(Sc.newdata,3)],0);
end

D2=timeonset(D2,Sc.time(1));

if(length(Sc.time)>1)
    D2 = fsample(D2, 1/diff(Sc.time(1:2)));
end

%%%%%%%%%%%%%%%
%% setup channel_labels, which will be a list of the relevant
% channels/sensors, for use later in the Fieldtrip leadfield calculations
       
if ~isempty(modalities) && strcmp(modalities{1},'EEG')   % changed by DM
    modality_meeg='EEG';
else
    modality_meeg='MEG';
end
    
if isempty(modalities)    
    chanind = 1:size(Sc.newdata,1);
else
    chanind = strmatch(modality_meeg, D.chantype);
end

if(size(Sc.newdata,4)>1)
    % TF
    D2 = transformtype(D2, 'TF');
    D2 = D2.frequencies(:, Sc.frequencies);
    D2(chanind,:,:,:)=permute(Sc.newdata,[1 4 2 3]);

else
    D2(chanind,:,:)=Sc.newdata;    
end

for coni=1:length(Sc.cond_list),
    D2=conditions(D2,coni,Sc.cond_list{coni});
end

if size(Sc.D,1) == size(Sc.newdata,1)
    % We aren't changing the number of channels, keep all information
    if(length(Sc.D.badchannels)>0)
        D2 = badchannels(D2, D.badchannels, ones(length(D.badchannels),1));
    end
    D2=D2.chantype(chanind, D.chantype);
    D2=D2.chanlabels(chanind, D.chanlabels);
elseif ~isempty(modalities)
    % We are using a subset of modalities, keep only their information
    D2=D2.chantype(chanind, D.chantype(chanind));
    D2=D2.chanlabels(chanind, D.chanlabels(chanind));
else
    % We are completely changing the channels, change type to LFP
    D2=D2.chantype(chanind, chantype);
end

D2.save;


