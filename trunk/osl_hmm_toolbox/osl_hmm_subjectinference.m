function hmm = osl_hmm_subjectinference(S)

if nargin < 2
    options = [];
end

if ~isfield(S,'D') || ~isa(spm_eeg_load(S.D),'meeg')
    error('S.D should contain a valid SPM object!');
end
    
    
options.mode        = ft_getopt(S,'mode','raw');
options.pcadim      = ft_getopt(S,'pcadim',40);
options.whiten      = ft_getopt(S,'whiten',1);
options.normalise   = ft_getopt(S,'normalise',1);
options.K           = ft_getopt(S,'K',8);
options.nreps       = ft_getopt(S,'nreps',1);


trl = 1; % can make this compatible with trialwise data later...


switch options.mode
    
    case 'raw'
        
        
        
        C = osl_cov(D(:,:,trl));
        
    case {'envelope','ds_envelope'}
        D = osl_hilbenv(struct('D',D,'winsize',0));
             
end

[allsvd,MixingMatrix] = eigdec(C,options.pcadim);
clear C

% Whitening:
if options.whiten
    MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
end

% Get PCA components:
hmmdata = zeros(size(MixingMatrix,1),Denv.nsamples);
blks = osl_memblocks(Denv,2);
for i = 1:size(blks,1)
    hmmdata(:,blks(i,1):blks(i,2)) = MixingMatrix * Denv(:,blks(i,1):blks(i,2));
end

% Infer HMM:
hmm = osl_hmm_infer(hmmdata,struct('K',options.K,'order',0,'Ninits',options.nreps,'Hz',D.fsample));
hmm.MixingMatrix = MixingMatrix;
        
        
end






switch options.mode
    
    case 'raw'
        
        current_montage = montage(D,'getindex');
        mont = montage(D,'getmontage',current_montage);
        
        D = montage(D,'switch',0);
        
        trl = 1;
        
        % Sensor space covariance:
        C = osl_cov(D(:,:,trl));
        
        % Apply weights:
        C = mont.tra*C*mont.tra';
        
        % Eigendecomposition:
        [allsvd,MixingMatrix] = eigdec(C,options.pcadim);
        clear C
        
        % Whitening:
        if options.whiten
            MixingMatrix = MixingMatrix * diag(1./sqrt(allsvd));
        end
        
        if options.normalise
            MixingMatrix = diag(sqrt(diag(C))) * MixingMatrix;
        end
        
        MixingMatrix = MixingMatrix';
        
        % Apply beamformer and eigenvectors in one shot:
        hmmdata =  (MixingMatrix * mont.tra) * D(:,:,:);
        
        % Infer HMM:
        hmm = osl_hmm_infer(hmmdata,struct('K',options.K,'order',0,'Ninits',options.nreps,'Hz',D.fsample));
        hmm.MixingMatrix = MixingMatrix;
        
    case {'envelope','ds_envelope'}
        Denv = osl_hilbenv(struct('D',D,'winsize',0));
        
        % Eigendecomposition:
        C = osl_cov(Denv);
        [allsvd,MixingMatrix] = eigdec(C,options.pcadim);
        clear C
        
        % Whitening:
        if options.whiten
            MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
        end
        
        % Get PCA components:
        hmmdata = zeros(size(MixingMatrix,1),Denv.nsamples);
        blks = osl_memblocks(Denv,2);
        for i = 1:size(blks,1)
            hmmdata(:,blks(i,1):blks(i,2)) = MixingMatrix * Denv(:,blks(i,1):blks(i,2));
        end
        
        % Infer HMM:
        hmm = osl_hmm_infer(hmmdata,struct('K',options.K,'order',0,'Ninits',options.nreps,'Hz',D.fsample));
        hmm.MixingMatrix = MixingMatrix;
               
end
        
end



%maps = osl_hmm_statemaps(hmm,hmmdata,1,'pcorr');
