function [ D, pcadims, pcadim_all, norm_vec, normalisation_used, fig_handles fig_names ] = normalise_sensor_data( S )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalise modalities using smallest eigenvalues or mean of eigs (i.e. the overall variance)
%% calculated using good channels and good trials, and over all woi


if strcmp(S.modalities{1},'EEG')
    modality_meeg='EEG';
else
    modality_meeg='MEGANY';
end

use_fixed_scaling=0;

D=S.D;

do_plots = S.do_plots;

badind = indchantype(D,modality_meeg,'BAD');

samples2use = S.samples2use;
trials = S.trials;

if S.pca_dim == 0
    S.pca_dim = numel(setdiff(find(strncmpi(D.chantype,modality_meeg,3)), badind));
end

%%%%%
% calc normalisation using noise variance
norm_vec=ones(length(D.chanlabels),1);

clear chanind;

for ff=1:length(S.modalities)
    
    % get good channels
    chanind{ff}= D.indchantype(S.modalities{ff},'good');
    
    if isempty(chanind{ff})
        error(['No good ' S.modalities{ff} ' channels were found.']);
    end    
end

for ff=1:length(S.modalities)
    
    % calc normalisation
    tmpdat=D(chanind{ff},find(samples2use),trials);
    %remove epoch means
    tmpdat=permute(tmpdat,[1 3 2]);
    tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
    tmpdat=permute(tmpdat,[1 3 2]);
    
    %tmpdat=D(chanind{ff},find(samples2use),trials);
    dat=reshape(tmpdat,length(chanind{ff}),(sum(samples2use)*length(trials)))';
    
    vs{ff}=var(dat);
    
    [pcadim, allsvd]=establish_dim(dat,S);
    
    % normalise based on method
    if ~isfield(S, 'normalise_method') || isempty(S.normalise_method),
        warning([mfilename ':NormaliseMethodNotSet'], ...
            'Normalisation method not set. Using min_eig. \n');
        S.normalise_method = 'min_eig';
    end%if
    
    switch lower(S.normalise_method)
        case 'min_eig'
            normalisation(ff)=sqrt(mean(allsvd(pcadim-5:pcadim)));
        case 'mean_eig'
            normalisation(ff)=sqrt(mean(allsvd(1:end)));
        case 'none'
            normalisation(ff)=1;
        otherwise
            error([mfilename ':NormaliseMethodNotRecognised'], ...
                'Normalisation method not recognised. \n');
    end;
    
    pcadims(ff)=pcadim;
    allsvds{ff}=allsvd;
    
    disp(['Dimensionality for modality ' S.modalities{ff} ' is: ' num2str(pcadim)]);
    disp(['Modality ' S.modalities{ff} ' has smallest sqrt(eig) = ' num2str(normalisation(ff))]);
    
end;

%%%%%
% look at eigenspectrum over all good channels and trials
chanindall = setdiff(find(strncmpi(D.chantype,modality_meeg,3)), badind);

tmpdat=D(chanindall,find(samples2use),trials);
%remove epoch means
tmpdat=permute(tmpdat,[1 3 2]);
tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
tmpdat=permute(tmpdat,[1 3 2]);

dat=reshape(tmpdat,length(chanindall),(sum(samples2use)*length(trials)))';

if S.force_pca_dim,% GC 2013-10-22
    pcaDim = S.pca_dim;
else
    pcaDim = rank(dat'*dat)-1;
end%if

[Apca,~,allsvd] = pca(dat,'numcomponents',pcaDim);
allsvd = allsvd(1:pcaDim);

fig_handles(1)=sfigure;
if ~do_plots
    set(fig_handles(1),'visible','off');
end
fig_names{1}='prenormalised_eigs';
subplot(2,2,1);plot(log(allsvd));title('Pre-normalised log eigenspectrum');

%%%%%
% look at channel variances
cols={'r','g','b'};

Apca_mods=[];
clear Apca_mods2;
for ff=1:length(S.modalities),
    
    % calc normalisation
    tmpdat=D(chanind{ff},find(samples2use),trials);
    
    %remove epoch means
    tmpdat=permute(tmpdat,[1 3 2]);
    tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
    tmpdat=permute(tmpdat,[1 3 2]);
    
    %tmpdat=D(chanind{ff},find(samples2use),trials);
    dat=reshape(tmpdat,length(chanind{ff}),(sum(samples2use)*length(trials)))';
    
    vs{ff}=var(dat);
    
    [c,ia,ib] = intersect(chanind{ff},chanindall);
    
    Apca_mods=[Apca_mods; Apca(ib,:)];
    
    Apca_mods2{ff}=Apca(ib,:);
    
    subplot(2,2,2);plot(log(allsvds{ff}),cols{ff});ho;
    
end;

for ff=1:length(S.modalities),
    subplot(2,2,4);plot(std(Apca_mods2{ff})./std(Apca_mods),cols{ff});ho;
end;

subplot(2,2,2);legend(S.modalities);title('Pre-normalised log eigenspectrum');

subplot(2,2,3);plot(spm_vec(vs));xlabel('Chan');title('Pre-normalised variances');
subplot(2,2,4);xlabel('PC');title('Pre-normalised std ratios');legend(S.modalities);


if(use_fixed_scaling),
    if sum(normalisation),
        normalisation=1e-13*normalisation/sum(normalisation);
    else % prevent dividing by zero
        % do nothing - sum is zero
    end%if
end;

%normalisation=ones(size(normalisation))*1e-13;

%% apply normalisation
for ff=1:length(S.modalities),
    
    if normalisation(ff),% non-zero
        norm_vec(chanind{ff})=ones(length(chanind{ff}),1)./normalisation(ff);
    else % prevent dividing by zero
        norm_vec(chanind{ff})=ones(length(chanind{ff}),1);
    end%if
    
    disp(['Modality ' S.modalities{ff} ' has data normalisation ' num2str(normalisation(ff))]);
    
end

normalisation_used=normalisation;

tra = zeros(length(indchantype(D,S.modalities)),D.nchannels);
tra(:,indchantype(D,S.modalities)) = diag(norm_vec(indchantype(D,S.modalities)));
D = add_montage(D,tra,'normalised_sensors',D.chanlabels(indchantype(D,S.modalities)));


%% establish dim of ALL normalised data

badind = indchantype(D,modality_meeg,'BAD');

chanindall = setdiff(find(strncmpi(D.chantype,modality_meeg,3)), badind);

% recalc chaninds
clear chanind;

for ff=1:length(S.modalities),
    
    % get good channels
    chanind{ff} = strmatch(S.modalities{ff}, D.chantype);
    chanind{ff} = setdiff(chanind{ff}, D.badchannels);
    if isempty(chanind{ff})
        error(['No good ' S.modalities{ff} ' channels were found.']);
    end    
end;

%dat=reshape(Dnew(chanindall,find(samples2use),trials),length(chanindall),(numel(find(samples2use))*length(trials)))';
tmpdat=D(chanindall,find(samples2use),trials);

%remove epoch means
tmpdat=permute(tmpdat,[1 3 2]);
tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
tmpdat=permute(tmpdat,[1 3 2]);
dat=reshape(tmpdat,length(chanindall),(sum(samples2use)*length(trials)))';

[pcadim_all allsvd_new]=establish_dim(dat,S);

% if(pcadim_all>=sum(pcadims))
%     pcadim_all=sum(pcadims)-1;
% end;

if(pcadim_all>=min(pcadims))
    pcadim_all=min(pcadims)-1;
end;

disp(['Overall dimensionality is: ' num2str(pcadim_all)]);


for ff=1:length(S.modalities),
    
    % calc normalisation
    tmpdat=D(chanind{ff},find(samples2use),trials);
    %remove epoch means
    tmpdat=permute(tmpdat,[1 3 2]);
    tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
    tmpdat=permute(tmpdat,[1 3 2]);
    %dat=reshape(Dnew(chanind{ff},find(samples2use),trials),length(chanind{ff}),(numel(find(samples2use))*length(trials)))';
    dat=reshape(tmpdat,length(chanind{ff}),(sum(samples2use)*length(trials)))';
    
    vs{ff}=var(dat);
    
    [pcadim allsvd]=establish_dim(dat,S);
    
    %normalisation(ff)=sqrt(mean(allsvd(pcadim-5:pcadim)))
    normalisation(ff)=sqrt(mean(allsvd(1:end)));
    
    pcadims(ff)=pcadim;
    allsvds{ff}=allsvd;
    
    disp(['Dimensionality for modality ' S.modalities{ff} ' is: ' num2str(pcadim)]);
    disp(['Modality ' S.modalities{ff} ' has smallest sqrt(eig) = ' num2str(normalisation(ff))]);
    
end;

%%%%%
% look at eigenspectrum over all good channels and trials
chanindall = setdiff(find(strncmpi(D.chantype,modality_meeg,3)), badind);

tmpdat=D(chanindall,find(samples2use),trials);
%remove epoch means
tmpdat=permute(tmpdat,[1 3 2]);
tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
tmpdat=permute(tmpdat,[1 3 2]);

dat=reshape(tmpdat,length(chanindall),(sum(samples2use)*length(trials)))';

if S.force_pca_dim,% GC 2013-10-22
    pcaDim = S.pca_dim;
else
    pcaDim = rank(dat'*dat)-1;
end%if

[Apca,~,allsvd] = pca(dat,'numcomponents',pcaDim);
allsvd = allsvd(1:pcaDim);

fig_handles(2)=sfigure;
if ~do_plots
    set(fig_handles(2),'visible','off');
end
fig_names{2}='normalised_eigs';
subplot(2,2,1);plot(log(allsvd));title('normalised log eigenspectrum');

%%%%%
% look at channel variances
cols={'r','g','b'};

Apca_mods=[];
clear Apca_mods2;
for ff=1:length(S.modalities),
    
    tmpdat=D(chanind{ff},find(samples2use),trials);
    %remove epoch means
    tmpdat=permute(tmpdat,[1 3 2]);
    tmpdat=tmpdat-repmat(mean(tmpdat,3),[1,1,sum(samples2use)]);
    tmpdat=permute(tmpdat,[1 3 2]);
    
    dat=reshape(tmpdat,length(chanind{ff}),(sum(samples2use)*length(trials)))';
    
    vs{ff}=var(dat);
    
    [c,ia,ib] = intersect(chanind{ff},chanindall);
    
    Apca_mods=[Apca_mods; Apca(ib,:)];
    
    Apca_mods2{ff}=Apca(ib,:);
    
    subplot(2,2,2);plot(log(allsvds{ff}),cols{ff});ho;
    
end;

for ff=1:length(S.modalities),
    subplot(2,2,4);plot(std(Apca_mods2{ff})./std(Apca_mods),cols{ff});ho;
end;

subplot(2,2,2);legend(S.modalities);title('normalised log eigenspectrum');

subplot(2,2,3);plot(spm_vec(vs));xlabel('Chan');title('normalised variances');
subplot(2,2,4);xlabel('PC');title('normalised std ratios');legend(S.modalities);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pcadim, allsvd]=establish_dim(dat,S)

% check for empty input
if isempty(dat)
    pcadim = 0;
    allsvd = [];
    return;
end

% setup pca rank
pcadim=S.pca_dim;

if S.force_pca_dim
    disp('Forcing PCA rank to be user-specified value.');
else
    pcadim_adapt = spm_pca_order(dat)-1;
    if((pcadim==-1 || pcadim>pcadim_adapt) && pcadim_adapt>1)
        pcadim = pcadim_adapt;
    end
    
end

[Apca,~,allsvd] = pca(dat,'numcomponents',pcadim);
allsvd = allsvd(1:pcadim);
min_eig2use = osl_check_eigenspectrum(allsvd, pcadim, 0);

if S.force_pca_dim
    if(min_eig2use<S.pca_dim)
        disp(['min_eig2use=' num2str(min_eig2use)]);
        disp(['S.pca_dim=' num2str(S.pca_dim)]);
        
        error('Dimensionality of data is less than the pca_dim being forced');
    end
else
    pcadim=min_eig2use;
end
