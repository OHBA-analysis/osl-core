function [ res ] = teh_spectral_nnmf( S )

% [ res ] = teh_spectral_nnmf( S )
%

try tmp=S.do_plots; catch, S.do_plots=1; end

res=[];

num_nodes=size(S.psds,4);
nfreqbins=size(S.psds,3);
nsubjects=size(S.psds,1);
NK=size(S.psds,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute coherency and auto spectra

coh_comps = zeros(size(S.psds,1),size(S.psds,2),nfreqbins,num_nodes,num_nodes); 
auto_spectra_comps = zeros(size(S.psds,1),size(S.psds,2),nfreqbins,num_nodes); 

for ss = 1:nsubjects
    for kk = 1:NK
        psd=squeeze(S.psds(ss,kk,:,:,:));
        for j=1:num_nodes,
            auto_spectra_comps(ss,kk,:,j) = psd(:,j,j);
            for l=1:num_nodes,                
                cjl = psd(:,j,l)./sqrt(psd(:,j,j) .* psd(:,l,l));
                coh_comps(ss,kk,:,j,l) = cjl;
            end
        end
    end
end

if 0
    % plot sanity check
    ii=37,jj=38;
    figure;
    subplot(221);plot(squeeze(abs(auto_spectra_comps(1,:,:,ii)))');ho;
    subplot(223);plot(squeeze(abs(coh_comps(1,:,:,ii,jj)))');ho;
    subplot(224);plot(squeeze(abs(auto_spectra_comps(1,:,:,jj)))');ho;
    subplot(221);plot(squeeze(abs(auto_spectra_comps(1,9,:,ii)))', 'LineWidth',3);ho;
    subplot(223);plot(squeeze(abs(coh_comps(1,9,:,ii,jj)))', 'LineWidth',3);ho;
    subplot(224);plot(squeeze(abs(auto_spectra_comps(1,9,:,jj)))', 'LineWidth',3);ho;
    legend(num2str([1:NK]'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NNMF ON PSD concat over states

if isfield(S,'nnmf_psd_specs') && isfield(S,'nnmf_psd_maps')
    res.nnmf_psd_specs=S.nnmf_psd_specs;
    res.nnmf_psd_maps=S.nnmf_psd_maps;
    
else    
    maxP=S.maxP;
    
    % concat over states    
    psdtmp=[];
    for kk=1:NK
        psdtmp=cat(2,psdtmp, squeeze(mean(abs(auto_spectra_comps(:,kk,:,:)),1)));    
    end

    if ~isfield(S,'fixed_psd_specs')
        [a b]=nnmf_mww(psdtmp,maxP,'replicates',500,'algorithm','als');
    else
        opt = statset('maxiter',1);    
        [a b]=nnmf_mww(psdtmp,size(S.fixed_psd_specs,1),'algorithm','als','w0',S.fixed_psd_specs','options',opt);
    end
    
    maps=[];
    specs=[];
    for pp=1:maxP
        ind=1;
        for kk=1:NK
            maps(kk,pp,:)=b(pp,ind:ind+num_nodes-1);
            ind=ind+num_nodes;
        end
        specs(pp,:)=a(:,pp)';
    end

    % put modes in order of increasing peak frequency
    [x,ind]=max(specs,[],2);
    [x,neworder_auto]=sort(ind);

    res.nnmf_psd_specs=specs(neworder_auto,:);
    res.nnmf_psd_maps=maps(:,neworder_auto,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Now do coh

clear a statset
for pp=1:S.maxP
    a(:,pp)=res.nnmf_psd_specs(pp,:)';
end

% concat coh over classes
psdtmp=[];
for kk=1:NK
    
    % extract lower diagonal only
    cohabs=squeeze(mean(abs(coh_comps(:,kk,:,:,:)),1));
    clear psd2;
    for pp=1:size(cohabs,1)
        psdf=triu(squeeze(cohabs(pp,:,:)),1);
        [inds]=(psdf~=0);
        psd2(pp,:)=psdf(inds);
    end
    num_per_kk=size(psd2,2);
    
    psdtmp=cat(2,psdtmp, psd2);
        
end


maxPcoh=S.maxPcoh;

if isfield(S,'nnmf_coh_specs') && isfield(S,'nnmf_coh_maps')
    res.nnmf_coh_specs=S.nnmf_coh_specs;
    res.nnmf_coh_maps=S.nnmf_coh_maps;
    
else   
    if ~isfield(S,'fixed_coh_specs')
        [anew bnew]=nnmf_mww(psdtmp,maxPcoh,'replicates',500,'algorithm','als');
    else
        opt = statset('maxiter',1);    
        [anew bnew]=nnmf_mww(psdtmp,size(S.fixed_coh_specs,1),'algorithm','als','w0',S.fixed_coh_specs','options',opt);
    end
         
    for pp=1:S.maxPcoh   
        coh_specs(pp,:)=anew(:,pp)';   
    end

    % put modes in order of increasing peak frequency
    [x,ind]=max(coh_specs,[],2);
    [x,neworder_auto]=sort(ind);

    res.nnmf_coh_specs=coh_specs(neworder_auto,:);

    %%%%
    % extract spatial info
    bnew_reordered=bnew(neworder_auto,:)+1e-10;
    
    res.nnmf_coh_maps=zeros(NK,S.maxPcoh,num_nodes,num_nodes);
    for pp=1:S.maxPcoh
        ind=1;
        for kk=1:NK        
            graph=bnew_reordered(pp,ind:ind+num_per_kk-1);
            graphmat=zeros(num_nodes,num_nodes);

            graphmat(inds)=graph;
            graphmat=(graphmat+graphmat')/2;

            res.nnmf_coh_maps(kk,pp,:,:)=graphmat;
            ind=ind+num_per_kk;
        end

    end
end


%%%%%%%%%%%%%%%%%%
%% plot specs
if S.do_plots
    figure;
    for pp=1:S.maxP   
        subplot(121);plot(res.nnmf_psd_specs(pp,:),get_cols(pp),'Linewidth',2);ho;
    end
    for pp=1:S.maxPcoh   
        subplot(122);plot(res.nnmf_coh_specs(pp,:),get_cols(pp),'Linewidth',2);ho;
    end
    legend(get_cols)
end

end

