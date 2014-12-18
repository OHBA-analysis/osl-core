function [hmm,maps] = osl_hmm_subjectinference(D,options)

if nargin < 2
    options = [];
end

options.mode    = ft_getopt(options,'mode','raw');
options.pcadim  = ft_getopt(options,'pcadim',40);
options.whiten  = ft_getopt(options,'whiten',1);
options.K       = ft_getopt(options,'K',8);
options.nreps   = ft_getopt(options,'nreps',1);



switch options.mode
    
    case 'raw'
        
        current_montage = montage(D,'getindex');
        mont = montage(D,'getmontage',current_montage);
        
        D = montage(D,'switch',0);
        
        trl = 1;
        
        % Sensor space covariance:
        C = cov(D(:,:,trl)');
        
        % Apply weights:
        C = mont.tra*C*mont.tra';
        
        % Eigendecomposition:
        [allsvd,MixingMatrix] = eigdec(C,options.pcadim);
        clear C
        
        % Whitening:
        if options.whiten
            MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
        end
        
        % Apply beamformer and eigenvectors in one shot:
        hmmdata =  (MixingMatrix * mont.tra) * D(:,:,:);
        
        
        % Infer HMM:
        hmm = osl_hmm_infer(hmmdata,struct('K',options.K,'order',0,'Ninits',options.nreps,'Hz',D.fsample));
        hmm.MixingMatrix = MixingMatrix;
        
    case {'envelope','ds_envelope'}
        Denv = osl_hilbenv(struct('D',D,'winsize',0));
        C = osl_cov(Denv);
% S.winsize
                
        
end
        
end



%maps = osl_hmm_statemaps(hmm,hmmdata,1,'pcorr');
