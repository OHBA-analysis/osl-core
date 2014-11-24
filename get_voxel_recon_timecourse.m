function [ dat wnorms wnorm wnorms_nai wnorm_nai ] = get_voxel_recon_timecourse( S )

% [ dat wnorms wnorm wnorms_nai wnorm_nai ] = get_voxel_recon_timecourse( S )
%
% MWW 2012

source_recon_results=S.source_recon_results;
chanind=S.chanind;
class_samples_inds=S.class_samples_inds;
voxind=S.voxind; %=first_level_results.mask_indices_in_source_recon(indind)

if ~strcmp(source_recon_results.recon_method,'none'), % not sensor space analysis
    if isfield(source_recon_results.BF.inverse.W,'MEG') % added by DM
        NK=length(source_recon_results.BF.inverse.W.MEG);
    else
        NK=length(source_recon_results.BF.inverse.W.EEG);
    end
else
    NK=1;
end;

Ntrials=size(class_samples_inds{1},3);
Ntpts=size(class_samples_inds{1},2);
Nfreqs=S.num_freqs;

dat=nan(Ntrials,Ntpts,Nfreqs);
wnorms=nan(Ntrials,Ntpts,Nfreqs);
wnorms_nai=nan(Ntrials,Ntpts,Nfreqs);
wnorm=zeros(NK,1);
wnorm_nai=zeros(NK,1);
weights=zeros(NK,length(chanind));

for kk=1:NK,

    if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis

        % for sensor space - just compute for all sensors
        weights(kk,voxind)=1;
    else
        
        if isfield(source_recon_results.BF.inverse.W,'MEG')  % added by DM
            weights(kk,:)=source_recon_results.BF.inverse.W.MEG{kk}{voxind};
        else
            weights(kk,:)=source_recon_results.BF.inverse.W.EEG{kk}{voxind};
        end
    end;

    if sum(isnan(squash(weights(kk,:))))==0,

        if ~strcmp(source_recon_results.recon_method,'none'),
            % this one represents the uncertainty and will be
            % applied equally to the data and regressors, is
            % therefore irrelevant if NK=1
            
            if isfield(source_recon_results.BF.inverse.W,'MEG')  % added by DM
                cov=source_recon_results.BF.features.C.MEG{kk};
            else
                cov=source_recon_results.BF.features.C.EEG{kk};
            end
            
            wnorm(kk)=weights(kk,:)*cov*weights(kk,:)';
            % this one represents a scaling of just the data
            wnorm_nai(kk)=weights(kk,:)*weights(kk,:)';   

        else
            wnorm(kk)=1;
            wnorm_nai(kk)=1;
        end;
    else,
        disp('NAN weights');
    end;
    
end;

for kk=1:NK,
    wkk=weights(kk,:);
    for tri=1:Ntrials, % indexes trials
        if ~isempty(S.timeinds{kk,tri}),
            dat(tri,S.timeinds{kk,tri},:)=reshape(wkk*S.sensor_data_sub{kk,tri},[1, length(S.timeinds{kk,tri}), Nfreqs]);

            wnorms(tri,S.timeinds{kk,tri},:)=wnorm(kk);                   
            wnorms_nai(tri,S.timeinds{kk,tri},:)=wnorm_nai(kk); 
        end;
    end;
end;

end

