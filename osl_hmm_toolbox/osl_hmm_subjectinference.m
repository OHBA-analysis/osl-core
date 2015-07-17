function hmm = osl_hmm_subjectinference(S)

if ~isfield(S,'D') || ~isa(spm_eeg_load(S.D),'meeg')
    error('S.D should contain a valid SPM object!');
end


S.mode        = ft_getopt(S,'mode','raw');
S.pcadim      = ft_getopt(S,'pcadim',40);
S.whiten      = ft_getopt(S,'whiten',1);
S.normalise   = ft_getopt(S,'normalise',1);

trl = 1; % can make this compatible with trialwise data later...

options.K      = ft_getopt(S,'K',8);
options.Ninits = ft_getopt(S,'nreps',1);


switch S.mode
    
    case 'raw'
        D = spm_eeg_load(S.D);
        options.zeromean = 1;
        
    case {'envelope','ds_envelope'}
        D = osl_hilbenv(struct('D',S.D,'winsize',0));
        options.zeromean = 0;
        
end

C = osl_cov(D);

% Normalisation:
if S.normalise
    % cov(Ax+a) = A*cov(x)*A'
    % xn = (x - mu)/sigma;
    % => cov(xn) = (1/sigma)*cov(x)*(1/sigma)';
    sigma = sqrt(osl_source_variance(D));
    C = diag(1./sigma)*C*diag(1./sigma);
end

[allsvd,MixingMatrix] = eigdec(C,S.pcadim);
MixingMatrix = MixingMatrix';

clear C

% Whitening:
if S.whiten
    MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix;
end

% Get PCA components:
hmmdata = zeros(size(MixingMatrix,1),D.nsamples);
blks = osl_memblocks(size(D),2);
for i = 1:size(blks,1)
    datblk = D(:,blks(i,1):blks(i,2),trl);
    if S.normalise
        datblk = normalise(datblk,2);
    end
    hmmdata(:,blks(i,1):blks(i,2)) = MixingMatrix * datblk;
end

% Infer HMM:
options.Hz      = D.fsample;
options.order   = 0;
options.initrep = 1;
options.zeromean = 1;
%options.initcyc = 10;
hmm = osl_hmm_infer(hmmdata,options);
hmm.MixingMatrix = MixingMatrix;


end

