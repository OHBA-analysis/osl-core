function [ dat wnorms wnorm wnorms_nai wnorm_nai weightsout ] = get_voxel_recon_timecourse_vector( S )

% [ dat wnorms wnorm wnorms_nai wnorm_nai ] = get_voxel_recon_timecourse_vector( S )
%
% Can handle scalar and vector recon
%
% MWW 2013

source_recon_results=S.source_recon_results;
chanind=S.chanind;
class_samples_inds=S.class_samples_inds;
voxind=S.voxind; %=first_level_results.mask_indices_in_source_recon(indind)
try, reduce_rank=S.reduce_rank; catch reduce_rank=1; end;

if ~strcmp(source_recon_results.recon_method,'none'), % not sensor space analysis
    if isfield(source_recon_results.BF.inverse,'MEG') % added by DM
        NK=length(source_recon_results.BF.inverse.MEG.class);
    else
        NK=length(source_recon_results.BF.inverse.EEG.class);
    end
else
    NK=1;
end;

Ntrials=size(class_samples_inds{1},3);
Ntpts=size(class_samples_inds{1},2);
Nfreqs=S.num_freqs;

% dipole rank
if ~strcmp(source_recon_results.recon_method,'none'), % sensor space analysis
    if isfield(source_recon_results.BF.inverse,'MEG') % added by DM
        Nrank=size(source_recon_results.BF.inverse.MEG.class{1}.W{1},1);    
    else
        Nrank=size(source_recon_results.BF.inverse.EEG.class{1}.W{1},1);    
    end
    if(reduce_rank>1)
        Nrank=reduce_rank;
    end;
    
else
    Nrank=1;
end;

dat=nan(Ntrials,Ntpts,Nfreqs,Nrank);
wnorms=nan(Ntrials,Ntpts,Nfreqs,Nrank);
wnorms_nai=nan(Ntrials,Ntpts,Nfreqs,Nrank);
wnorm=zeros(NK,Nrank);
wnorm_nai=zeros(NK,Nrank);
weights=zeros(NK,Nrank,length(chanind));


for kk=1:NK,
    if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis
        % for sensor space - just compute for all sensors
        weights(kk,1,voxind)=1;
    else
        if isfield(source_recon_results.BF.inverse,'MEG')  % added by DM
            w=source_recon_results.BF.inverse.MEG.class{kk}.W{voxind};
        else
            w=source_recon_results.BF.inverse.EEG.class{kk}.W{voxind};
        end
                
        if(reduce_rank>1)
            tmp=w*w';   
            [u, ~] = svd(tmp,'econ');                              
            w = u(:,1:reduce_rank)'*w;
        end;
        
        weights(kk,:,:)=w;
    end;
    
    if sum(isnan(squash(weights(kk,:,:))))==0,
        
        if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis
            cov=eye(size(weights,3));
        else
            if isfield(source_recon_results.BF.inverse,'MEG')  % added by DM
                cov=source_recon_results.BF.features.MEG.class{kk}.C;
            else
                cov=source_recon_results.BF.features.EEG.class{kk}.C;
            end
            
            if NK>1
                % compute covs over other states
                cov_mkk{kk}=zeros(size(cov));                
                for kk2=1:NK,
                    if kk~=kk2,
                        if isfield(source_recon_results.BF.inverse,'MEG')  % added by DM
                            covtmp=source_recon_results.BF.features.MEG.class{kk2}.C;
                        else
                            covtmp=source_recon_results.BF.features.EEG.class{kk2}.C;
                        end
                        cov_mkk{kk}=cov_mkk{kk}+covtmp;
                    end;
                end;
                cov_mkk{kk}=cov_mkk{kk}/(NK-1);
            else
                cov_mkk{kk}=eye(size(weights,3));
            end;
        end;
        
        
        weightskk=permute(weights(kk,:,:),[2 3 1]);
                
        wnorm(kk)=trace(weightskk*cov*weightskk');
        %wnorm(kk)=trace(weightskk*cov_mkk{kk}*weightskk');
        wnorm_nai(kk)=trace(weightskk*weightskk');

    else,
        disp('NAN weights');
    end;
    
end;

weightsout{kk}=[];

if(1),
    for kk=1:NK,

        weightskk=permute(weights(kk,:,:),[2 3 1]);            
        weightsout{kk}=weightskk;

        for tri=1:Ntrials, % indexes trials

            if ~isempty(S.timeinds{kk,tri}),
                tmp=weightskk*S.sensor_data_sub{kk,tri}; % nrank x (ntpts x nfreqs)
                tmp2=reshape(tmp,[Nrank,length(S.timeinds{kk,tri}),Nfreqs]); % nrank x ntpts x nfreqs
                %dat(tri,S.timeinds{kk,tri},:,:)=reshape(tmp,[1, length(S.timeinds{kk,tri}),Nfreqs,Nrank];);
                dat(tri,S.timeinds{kk,tri},:,:)=permute(tmp2,[4 2 3 1]); % Ntrials x Ntpts x Nfreqs x Nrank

                if(1)
                    wnorms(tri,S.timeinds{kk,tri},:,:)=wnorm(kk);
                    wnorms_nai(tri,S.timeinds{kk,tri},:,:)=wnorm_nai(kk);
                else
                    wnorms(tri,S.timeinds{kk,tri},:,:)=mean(wnorm);
                    wnorms_nai(tri,S.timeinds{kk,tri},:,:)=mean(wnorm_nai);
                end;
            end;
        end;

    end;
else
    weights_s=zeros(Ntpts,length(chanind),Nfreqs,Nrank);
    
    for tri=1:Ntrials, % indexes trials

        for kk=1:NK,
        
            weightskk=permute(weights(kk,:,:),[2 3 1]); % nrank x chans           

            if ~isempty(S.timeinds{kk,tri}),
                tmp=weightskk; % nrank x chans
                  
                tmp2=repmat(weightskk,[1,1,1,length(S.timeinds{kk,tri}),Nfreqs]);
                
                weights_s(S.timeinds{kk,tri},:,:,:)=permute(tmp2,[4 2 3 1]);
                
                wnorms(tri,S.timeinds{kk,tri},:,:)=mean(wnorm);
                wnorms_nai(tri,S.timeinds{kk,tri},:,:)=mean(wnorm_nai);

            end;
        end;
        
        for tt=1:Ntpts,
            dat(tri,tt,:,:)=weights_s(tt,:)*S.sensor_data_tf(:,tt,tri);
        end;
        
    end;
end;

end

