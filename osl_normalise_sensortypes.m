function D = osl_normalise_sensortypes(S)
% Normalise Elekta sensortypes using minimum eigenvalue normalisation
%
% D = osl_normalise_sensortypes(S)
%
% REQUIRED INPUTS:
%
% S.D        - SPM MEG object filename
% 
% OPTIONAL INPUTS:
%
% S.method   - normalisation method to use ['min,'mean'], default 'min'
%
% S.pcadim   - set the PCA dimensionality to use, default []
%
% S.prefix   - prefix to append to file name, default ''
%
% S.do_plots - produce diagnostic plots [0/1], default 0
%
% Adam Baker 2015

D = spm_eeg_load(S.D);

% Check SPM File Specification:
try
    S.D = char(S.D);
    [pathstr,filestr] = fileparts(S.D);
    S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
    D = spm_eeg_load(S.D);
catch
    error('SPM file specification not recognised or incorrect');
end

S.method   = ft_getopt(S,'method','min');
S.pcadim   = ft_getopt(S,'pcadim',[]);
S.prefix   = ft_getopt(S,'prefix','');
S.do_plots = ft_getopt(S,'do_plots',0);

modalities = {'MEGPLANAR','MEGMAG'};

if ~numel(intersect(D.chantype,modalities)) == 2
    error('SPM object must contain MEGPLANAR and MEGMAG modalities');
end

normalisation = ones(1,length(indchantype(D,modalities)));

if S.do_plots 
    figure
    linecolor = {'r','b'};
end

% Compute covariance of each sensor type
for modality = modalities
    
    % Good channels:
    channel_inds = indchantype(D,modality,'GOOD');
    
    % Select data:
    data = reshape(D(channel_inds,:,:),length(channel_inds),[]);
    
    % Good samples & trials:
    goodsamples = ~all(badsamples(D,channel_inds,':',':'));
    goodsamples = reshape(goodsamples,1,[]);
    
    % Remove bad samples:
    data = data(:,goodsamples);
    
    % Estimate normalisation
    C = osl_cov(data);
    
    if isempty(S.pcadim)
        rankEst = estimate_rank(C);
    else
        rankEst = S.pcadim;
    end
    
    eigenvectors = svd(C);
    
    switch S.method
        case 'mean'
            normalisation(channel_inds) = 1./sqrt(mean(eigenvectors(1:rankEst)));
        case 'min'
            normalisation(channel_inds) = 1./sqrt(eigenvectors(rankEst));
    end
    
    
    if S.do_plots
        subplot(1,2,1); title('Unnormalised eigenvectors'); hold on;
        plot(log(eigenvectors),'linewidth',2,'color',linecolor{strcmp(modality,modalities)})
        stem(rankEst,log(eigenvectors(rankEst)),'--','color',linecolor{strcmp(modality,modalities)},'marker','none')
        subplot(1,2,2); title('Normalised eigenvectors'); hold on;
        plot(log(eigenvectors.*normalisation(channel_inds).^2'),'linewidth',2,'color',linecolor{strcmp(modality,modalities)})
    end
    
end


% Apply normalisation via SPM montage:
montage             =  [];
montage.tra         =  diag(normalisation);
montage.labelnew    =  D.chanlabels(indchantype(D,modalities));
montage.labelorg    =  D.chanlabels(indchantype(D,modalities));

[~,indx] = ismember(montage.labelnew,D.sensors('MEG').label);

montage.chanunitnew =  D.sensors('MEG').chanunit(indx);
montage.chanunitorg =  D.sensors('MEG').chanunit(indx);
montage.chantypenew =  lower(D.sensors('MEG').chantype(indx));
montage.chantypeorg =  lower(D.sensors('MEG').chantype(indx));

Sm = [];
Sm.D              =  fullfile(D.path,D.fname);
Sm.montage        =  montage;
Sm.keepothers     =  true;
Sm.updatehistory  =  1;
Sm.prefix         = S.prefix;

if ~isempty(S.prefix)
    D = spm_eeg_montage(Sm);
else
    D = osl_spmfun(@spm_eeg_montage,Sm);
end
    


end


function rankEst = estimate_rank(C)
% Matlab's rank.m and SPM's spm_pca_order.m (Minka's approximation) seem to
% estimate the rank of Maxfiltered or AFRICAed data poorly... Please feel
% free to augment this function with a better estimation approach.

minRank = 40;

eigDiff = diff(log(svd(C)));

rankVec = eigDiff<10*median(eigDiff);
rankVec(1:minRank) = false;

rankEst = find(rankVec,1,'first');

if isempty(rankEst)
    rankEst = rank(C);
end

end
