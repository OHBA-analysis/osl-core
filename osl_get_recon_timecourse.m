function [ dat S ] = osl_get_recon_timecourse( S )

% [ dat ] = osl_get_recon_timecourse( S )
% [ dat S ] = osl_get_recon_timecourse( S )
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
% S.D_block.data - caches block of recon data for voxel indexes:
%                   S.D_block.from : S.D_block.from+S.D_block.size-1
%                   If not provided, or requested S.index is outside the 
%                   block, then will be reconstructed from S.D.
% S.D_block.from - Voxel index for start of block. If not provided, or 
%                   requested S.index is outside the block, then this will 
%                   be determined from S.index.
% S.D_block.size - No. of voxels in block. Defaults to 500
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

if ~isfield(S,'D_block')
    S.D_block=[];
end;
if ~isfield(S.D_block,'from')
    S.D_block.from=[];
end;
if ~isfield(S.D_block,'size')
    S.D_block.size=500;
end;  
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% establish dims
if ~isempty(S.D.nfrequencies),  
    nfreqs=S.D.nfrequencies;
else
    nfreqs=1;
end;

ntrials=S.D.ntrials;
ntpts=S.D.nsamples;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% does S.D contain a Class channel?
if ~isfield(S,'classchanind')
    classchanind=find(strcmp(S.D.chanlabels,'Class'));
    S.classchanind=classchanind;
else
    classchanind=S.classchanind;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if needed setup the S.class_timeinds
if ~isempty(classchanind),
    
    nclasses=max(squash(S.D(classchanind, :, :, 1)));

    if ~isfield(S,'class_timeinds'),

        %% establish the time indices for each class by using the 'Class' channel in S.D 
        S.class_timeinds=cell(nclasses,ntrials);

        for kk=1:nclasses,
            class_samples_inds_recon = (S.D(classchanind, :, :, 1)==kk);

            for tri=1:ntrials,            
                S.class_timeinds{kk,tri}=find(class_samples_inds_recon(1,:, tri)); % time indices for class kk
            end;
        end;        

    end;  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if S.index is outside the current block then create S.D_block from S.D

% is S.index outside the current block
use_current_block=0;
try
    use_current_block = (S.index >= S.D_block.from  && S.index <= S.D_block.from + S.D_block.size - 1);
catch
    use_current_block=0;
end;

if ~use_current_block

    % establish S.D_block.from to place block of size S.D_block.size around S.index
    S.D_block.from = floor((S.index-1)/S.D_block.size)*S.D_block.size+1;

    % do sanity check
    use_current_block = (S.index >= S.D_block.from  && S.index <= S.D_block.from + S.D_block.size - 1);
    if ~use_current_block, error('Something has gone wrong'); end;

    % construct new D_block from S.D 
    if ~isfield(S.D,'parcellation')
        D=montage(S.D,'switch',S.montage_index); 
    else
        D=S.D;
    end
    
    to=min(S.D_block.from+S.D_block.size-1,size(D,1));
    
    if isempty(classchanind),  
        %% there is no Class channel - so do a straightforward recon       
        S.D_block.data=D(S.D_block.from:to,:,:,:);
    else
        
        %% there is a Class channel               
        start=S.montage_index;

        S.D_block.data=nan(to-S.D_block.from+1,ntpts,ntrials,nfreqs);
        for kk=1:nclasses,
            % switch to recon montage for this class    
            if ~isfield(S.D,'parcellation')
                D=montage(S.D,'switch',start-1+kk);
            else
                D=S.D;
            end
            
            for ff=1:nfreqs,
                for tri=1:ntrials, 
                    if ~isempty(S.class_timeinds{kk,tri}),
                        S.D_block.data(:,S.class_timeinds{kk,tri},tri,ff)=D(S.D_block.from:to,S.class_timeinds{kk,tri},tri,ff); 
                    end;
                end;
            end;
        end;
        
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% get recon for voxel S.index from cached block of recon data

dat=S.D_block.data(S.index-S.D_block.from+1,:,:,:);



