function [ D2 ] = osl_concat_spm_eeg_chans( Sc )

% [ D2 ] = osl_concat_spm_eeg_chans( Sc )
%
% Clones a passed in Sc.D and concats passed in channels
% Needs:
% Sc.D

% Sc.newchandata (3D matrix of data (channels, samples, trials))
% Sc.newchanlabels (list of str names (channels))
% Sc.newchantype (list of str types (channels))

if(size(Sc.newchandata,4)>1)
    error('Only implemented for time domain data');
end;

if ~isfield(Sc, 'newname'),
    Sc.newname=[Sc.D.path '/concat' Sc.D.fname];
end

D2=clone(Sc.D,Sc.newname,[size(Sc.D,1)+size(Sc.newchandata,1),size(Sc.newchandata,2),size(Sc.newchandata,3)],0);

chanind=(size(Sc.D,1)+1):(size(Sc.D,1)+size(Sc.newchandata,1));

D2(1:size(Sc.D,1),:,:)=Sc.D(:,:,:);

D2(chanind,:,:)=Sc.newchandata;

D2 = chanlabels(D2, 1:size(Sc.D,1), Sc.D.chanlabels);
D2 = chantype(D2, 1:size(Sc.D,1), Sc.D.chantype);
D2 = units(D2, 1:size(Sc.D,1), units(Sc.D, 1:size(Sc.D,1)));

D2(chanind,:,:)=Sc.newchandata;
for ii=1:length(chanind),
    D2 = chanlabels(D2, chanind(ii), Sc.newchanlabels{ii});
    D2 = chantype(D2, chanind(ii), Sc.newchantype{ii});

end;

if(1)
    if ~isempty(Sc.D.badchannels)
        badchan = Sc.D.chanlabels(Sc.D.badchannels);
        badchanind = spm_match_str(D2.chanlabels, badchan);
        D2 = badchannels(D2, badchanind, 1);
    end
end;

D2.save;


