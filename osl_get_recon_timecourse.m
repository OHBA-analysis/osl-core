function [ dat class_timeinds ] = osl_get_recon_timecourse( S )

% [ dat ] = osl_get_recon_timecourse( S )
% [ dat class_timeinds ] = osl_get_recon_timecourse( S )
%  
% Reconstructs time courses using montages in S.D
%
% If S.D contains a 'Class' channel then this will combine recons from
% multiple classes.
%
% S.D - SPM MEEG object with montages (needs one for each class)
% 
% S.index - index of voxel to recon
%
% optional:
%
% S.class_timeinds - nclasses x ntrials cell matrix, where each cell contains a
% logical vector of length ntpts indicating membership for each class. 
% If not passed in then this will get created from the 'Class' channel in S.D.
% 
% S.montage_index - index of montage to use (if multi-class then this is
% the index of the the montage for the first class, i.e. the montage indexes 
% for the NK classes should correspond to the indexes 
% S.montage_index:S.montage_index+NK-1). Default value is 1.
%
% returns:
%
% dat - ntrials x ntpts x nfreq
%
% MWW 2014

S = ft_checkopt(S,'index',{'single','double'});

try
    S = ft_checkopt(S,'montage_index',{'single','double'});
catch 
    warning('montage_index, setting to 1')
    S = ft_setopt(S,'montage_index',1);
end

dat=[];
class_timeinds=[];

%% establish dims
if ~isempty(S.D.nfrequencies),  
    nfreqs=S.D.nfrequencies;
else
    nfreqs=1;
end;

ntrials=S.D.ntrials;
ntpts=S.D.nsamples;

%% setup container for output
dat=nan(ntrials,ntpts,nfreqs);

%% does S.D contain a Class channel?
classchanind=find(strcmp(S.D.chanlabels,'Class'));

if ~isempty(classchanind),
    
    %% there is no Class channel - so do a straightforward recon
    
    % switch to recon montage
    D=montage(S.D,'switch',S.montage_index);
    
    % recon
    dat=D(S.index,:,:,:);
    
else
    
    %% there is a Class channel
    nclasses=length(S.class_samples_inds);
    
    %% if needed setup the S.class_timeinds
    if ~isfield(S,'class_timeinds'),

        %% establish the time indices for each class by using the 'Class' channel in S.D 
        class_timeinds=cell(nclasses,ntrials);

        for kk=1:nclasses,
            class_samples_inds_recon = (S.D(classchanind, :, :, 1)==kk);

            for tri=1:ntrials,            
                class_timeinds{kk,tri}=find(class_samples_inds_recon(1,:, tri)); % time indices for class kk
            end;
        end;
    else
        class_timeinds=S.class_timeinds;
    end;    
    
    %% specify the index of the montage for the first class
    start=S.montage_index;

    %% recon
    for kk=1:nclasses,

        % switch to recon montage for this class    
        D=montage(S.D,'switch',start-1+kk);

        for ff=1:nfreqs,
            for tri=1:ntrials, 
                if ~isempty(class_timeinds{kk,tri}),
                    dat(tri,class_timeinds{kk,tri},ff)=D(S.index,class_timeinds{kk,tri},tri,ff); 
                end;
            end;
        end;

    end;

end;




