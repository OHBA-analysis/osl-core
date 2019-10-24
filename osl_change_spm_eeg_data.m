function [ D2 ] = osl_change_spm_eeg_data( Sc )

% [ D2 ] = osl_change_spm_eeg_data( Sc )
%
% Clones a passed in Sc.D and fills it with the passed in
% data. 
%
% Output D2 object will have number of chans equal to size(Sc.newdata,1).
%
% Number of channels in Sc.newdata needs to be the same size as chans of 
% type Sc.modality in Sc.D. Otherwise the channel info will be lost (and
% chantype will be set to Sc.chantype)
%
%
% Needs:
% Sc.D 
% Sc.newdata (3D matrix of data (channels, samples, trials))
% Sc.newname (filename - should end with suffix .dat)
%
% Optional:
% Sc.time (vector of time points)
% Sc.cond_list
% Sc.frequencies (vector of freqs)
% Sc.modalities ('MEG' or 'EEG'. Default leaves this empty and all chans in 
%                 D will be used)
% Sc.remove_montages (default is true)
% Sc.chantype (only gets used if Sc.newdata is not the same size as chans of 
%              type Sc.modality in Sc.D. Default is 'LFP')
% 
% MWW

try, tmp=Sc.cond_list; catch, Sc.cond_list=Sc.D.conditions; end
try, modalities=Sc.modalities; catch, modalities=[];  end
try, Sc.time=Sc.time; catch, Sc.time=Sc.D.time; end
try, Sc.remove_montages=Sc.remove_montages; catch, Sc.remove_montages=true; end
try, chantype=Sc.chantype; catch, chantype='LFP';  end

% remove montages
if Sc.remove_montages
    D=montage(Sc.D,'remove',1:montage(Sc.D,'getnumber'));
else
    D=Sc.D;
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
    % use all channels in D
    Dchanind = 1:size(D,1);
else
    Dchanind = strmatch(modality_meeg, D.chantype);
end

% Output D2 object will alwyas have number of chans equal to size(Sc.newdata,1).
D2chanind=1:size(Sc.newdata,1);

%%%%%%%%%%%%%%%
%% Number of channels in Sc.newdata needs to be the same size as chans of 
%% type Sc.modality in Sc.D. Otherwise all channel info will be lost
lose_channel_info=false;
if length(Dchanind) ~= size(Sc.newdata,1)
    warning('Number of channels in new data is NOT the same size as chans of type Sc.modality in Sc.D. All channel info will be lost, and passed in chantype will be used instead');
    lose_channel_info=true;
end

if(size(Sc.newdata,4)==1)
    %TF    
    D2=clone(D,Sc.newname,[size(Sc.newdata,1),size(Sc.newdata,2),size(Sc.newdata,3)],2);
else
    D2=clone(D,Sc.newname,[size(Sc.newdata,1),size(Sc.newdata,4),size(Sc.newdata,2),size(Sc.newdata,3)],2);
end

D2=timeonset(D2,Sc.time(1));

if(length(Sc.time)>1)
    D2 = fsample(D2, 1/diff(Sc.time(1:2)));
end

if(size(Sc.newdata,4)>1)
    % TF
    D2 = transformtype(D2, 'TF');
    D2 = D2.frequencies(:, Sc.frequencies);
    D2(D2chanind,:,:,:)=permute(Sc.newdata,[1 4 2 3]);

else
    D2(D2chanind,:,:)=Sc.newdata;    
end

for coni=1:length(Sc.cond_list),
    D2=conditions(D2,coni,Sc.cond_list{coni});
end

if lose_channel_info
    D2=D2.chantype(D2chanind, chantype);
else

    % keep channel information
    if(length(D.badchannels)>0)            
        D2 = badchannels(D2, find(D.badchannels(Dchanind)), ones(length(D.badchannels),1));
    end
    D2=D2.chantype(D2chanind, D.chantype(Dchanind));
    D2=D2.chanlabels(D2chanind, D.chanlabels(Dchanind));
    
end

if size(D2,3)==1
    D2=type(D2,'continuous');
end

D2.save;


