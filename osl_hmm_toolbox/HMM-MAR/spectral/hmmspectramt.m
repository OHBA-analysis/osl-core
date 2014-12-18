function fit = hmmspectramt(X,T,options)
%
% Computes nonparametric (multitaper) power, coherence, phase and PDC, with intervals of
% confidence. Also obtains mean and standard deviation of the phases.
% PDC is approximated: From MT cross-spectral density, it obtains the transfer-function
% using the Wilson-Burg algorithm and then computes PDC.
% Intervals of confidence are computed using the jackknife over trials and tapers
%
%
% INPUTS:
% X: the data matrix, with all trials concatenated
% T: length of each trial
% options: include the following fields
%   .Gamma: Estimated posterior probabilities of the states (default: all ones)
%   .tapers: A numeric vector [TW K] where TW is the
%       time-bandwidth product and K is the number of
%       tapers to be used (less than or equal to
%       2TW-1)
%   .win: number of time points per non-overlapping window
%   .Fs: Sampling frequency
%   .fpass: Frequency band to be used [fmin fmax] (default [0 fs/2])
%   .p: p-value for computing jackknife confidence intervals (default 0)
%   .numIterations: no. iterations for the Wilson algorithm (default: 100)
%   .tol: tolerance limit (default: 1e-18)
%
% OUTPUTS:
%
% fit is a list with K elements, each of which contains:
% fit.state(k).psd: cross-spectral density (Nf x N x N)
% fit.state(k).ipsd inverse power spectral density matrix [Nf x N x N]
% fit.state(k).coh: coherence (Nf x N x N)
% fit.state(k).pcoh partial coherence (Nf x N x N)
% fit.state(k).pdc: partial directed coherence estimate (Nf x N x N)
% fit.state(k).phase: phase of the coherence, in degrees (Nf x N x N)
% fit.state(k).psderr: interval of confidence for the cross-spectral density (2 x Nf x N x N)
% fit.state(k).coherr: interval of confidence for the coherence (2 x Nf x N x N)
% fit.state(k).pdcerr: interval of confidence for the partial directed coherence (2 x Nf x N x N)
% fit.state(k).sdphase: jackknife standard deviation of the phase
% fit.state(k).f: frequencies (Nf x 1)
%
% Author: Diego Vidaurre, OHBA, University of Oxford (2014)

ndim = size(X,2); 
if ~isfield(options,'Gamma'), Gamma = ones(sum(T),1); 
else Gamma = options.Gamma; end
if ~isfield(options,'p'), options.p = 0; end
if ~isfield(options,'removezeros'), options.removezeros = 0; end
if ~isfield(options,'rlowess'), options.rlowess = 0; end
if ~isfield(options,'numIterations'), options.numIterations = 100; end
if ~isfield(options,'tol'), options.tol = 1e-18; end
if ~isfield(options,'pad'), options.pad = 0; end;
if ~isfield(options,'Fs'), options.Fs=1; end;
if ~isfield(options,'fpass'),  options.fpass=[0 options.Fs/2]; end;
if ~isfield(options,'tapers'), options.tapers = [4 7]; end;
if ~isfield(options,'win'), win = min(T); else win = options.win; end
K = size(Gamma,2);

if options.p>0, options.err = [2 options.p]; end
Fs = options.Fs;
fit = {};

nfft = max(2^(nextpow2(win)+options.pad),win);
[f,findx]=getfgrid(Fs,nfft,options.fpass);
Nf = length(f);
tapers=dpsschk(options.tapers,win,Fs); % check tapers
ntapers = options.tapers(2);

% remove the exceeding part of X (with no attached Gamma)
order = (size(X,1) - size(Gamma,1)) / length(T);
X2 = [];
for in=1:length(T),
    t0 = sum(T(1:in-1));
    X2 = [X2; X(t0+1+order:t0+T(in),:)];
    T(in) = T(in) - order;
end
X = X2; clear X2; 

for k=1:K
    
    Xk = X .* repmat(Gamma(:,k),1,ndim);
    if options.removezeros
        nonzero = Gamma~=0;
        Xk = Xk(nonzero,:);
        Tk = T;
        for in=1:length(T),
            t0 = sum(T(1:in-1));
            Tk(in) = sum(Gamma(t0+1:t0+T(in))~=0);
        end
    else
        Tk = T;
    end
    
    % Multitaper Cross-frequency matrix calculation
    psdc = zeros(Nf,ndim,ndim,length(Tk)*ntapers);
    for in=1:length(Tk)
        Nwins=round(Tk(in)/win); % that means that a piece of data is going to be included as a window only if it's long enough
        t0 = sum(Tk(1:in-1));
        Xki = Xk(t0+1:t0+Tk(in),:);
        for iwin=1:Nwins,
            ranget = (iwin-1)*win+1:iwin*win;
            if ranget(end)>Tk(in),
                nzeros = ranget(end) - Tk(in);
                nzeros2 = floor(nzeros/2) * ones(2,1);
                if sum(nzeros2)<nzeros, nzeros2(1) = nzeros2(1)+1; end
                ranget = ranget(1):Tk(in);
                Xwin=[zeros(nzeros2(1),ndim); Xki(ranget,:); zeros(nzeros2(2),ndim)]; % padding with zeroes
            else
                Xwin=Xki(ranget,:);
            end
            J=mtfftc(Xwin,tapers,nfft,Fs); % use detrend on X?
            for tp=1:ntapers,
                Jik=J(findx,tp,:);
                for j=1:ndim,
                    for l=1:ndim,
                        psdc(:,j,l,(in-1)*ntapers+tp) = psdc(:,j,l,(in-1)*ntapers+tp) + mean(conj(Jik(:,:,j)).*Jik(:,:,l),2) / Nwins;
                    end
                end
            end
        end
    end
    
    % coherence
    psd = mean(psdc,4); ipsd = zeros(Nf,ndim,ndim);
    for ff=1:Nf, ipsd(ff,:,:) = inv(permute(psd(ff,:,:),[3 2 1])); end
    coh = zeros(Nf,ndim,ndim); phase = zeros(Nf,ndim,ndim);
    for j=1:ndim, 
        for l=1:ndim,
            cjl = psd(:,j,l)./sqrt(psd(:,j,j) .* psd(:,l,l));
            coh(:,j,l) = cjl;
            pcoh(:,j,l) = -ipsd(:,j,l)./sqrt(ipsd(:,j,j) .* ipsd(:,l,l));
            phase(:,j,l) = angle(cjl); %atan(imag(rkj)./real(rkj));
        end
    end
    
    pdc = subrutpdc(psd,options.numIterations,options.tol);
    psderr = []; coherr = []; pdcerr = [];
    
    if options.p>0 % jackknife
        if length(Tk)*options.tapers(2)<5,  error('You need at least 5 trials*tapers to compute error bars\n'); end
        jksamples = size(psdc,4);
        SA = []; cohA = []; pdcA = []; dcA = [];
        atanhcohA = [];
        atanhpdcA = [];
        phasefactorA = [];
        for in=1:jksamples;
            %fprintf('%d of %d \n',in,size(psdc,4))
            indxk=setdiff(1:jksamples,in);
            Si = mean(psdc(:,:,:,indxk),4);
            cohi = zeros(Nf,ndim,ndim);
            atanhcohi = zeros(Nf,ndim,ndim);
            atanhpdci = zeros(Nf,ndim,ndim);
            phasefactori = zeros(Nf,ndim,ndim);
            pdci = subrutpdc(Si,options.numIterations,options.tol);
            for j=1:ndim,
                for l=1:ndim,
                    cjl = Si(:,j,l)./sqrt(Si(:,j,j) .* Si(:,l,l));
                    cohi(:,j,l) = abs(cjl);
                    atanhcohi(:,j,l) = sqrt(2*jksamples-2)*atanh(cohi(:,j,l));
                    atanhpdci(:,j,l) = sqrt(2*jksamples-2)*atanh(pdci(:,j,l));
                    phasefactori(:,j,l)=cjl./cohi(:,j,l);
                end
            end
            SA = cat(4,SA,Si); cohA = cat(4,cohA,cohi);
            atanhcohA = cat(4,atanhcohA,atanhcohi);
            atanhpdcA = cat(4,atanhpdcA,atanhpdci);
            phasefactorA = cat(4,phasefactorA,phasefactori);
            pdcA = cat(4,pdcA,pdci); dcA = cat(4,pdcA,pdci);
        end
        dof = jksamples-1;
        tcrit=tinv(1-options.p/2,dof); % Inverse of Student's T cumulative distribution function
        % psd errors
        sigmaS = tcrit * sqrt(dof)*std(log(SA),1,4);
        psderr = zeros([2 size(psd)]);
        psderr(1,:,:,:) = psd .* exp(-sigmaS);
        psderr(2,:,:,:) = psd .* exp(sigmaS);
        % coh errors
        atanhcoh=sqrt(2*jksamples-2)*atanh(coh); % z
        sigma12=sqrt(jksamples-1)*std(atanhcohA,1,4); % Jackknife estimate std(z)=sqrt(dim-1)*std of 1-drop estimates
        Cu=atanhcoh+tcrit*sigma12;
        Cl=atanhcoh-tcrit*sigma12;
        coherr = zeros([2 size(coh)]);
        coherr(1,:,:,:) = max(tanh(Cl/sqrt(2*jksamples-2)),0); % This ensures that the lower confidence band remains positive
        coherr(2,:,:,:) = tanh(Cu/sqrt(2*jksamples-2));
        % pdc errors
        atanhpdc=sqrt(2*jksamples-2)*atanh(pdc); % z
        sigma12=sqrt(jksamples-1)*std(atanhpdcA,1,4); % Jackknife estimate std(z)=sqrt(dim-1)*std of 1-drop estimates
        Cu=atanhpdc+tcrit*sigma12;
        Cl=atanhpdc-tcrit*sigma12;
        pdcerr = zeros([2 size(pdc)]);
        pdcerr(1,:,:,:) = max(tanh(Cl/sqrt(2*jksamples-2)),0); % This ensures that the lower confidence band remains positive
        pdcerr(2,:,:,:) = tanh(Cu/sqrt(2*jksamples-2));
        % std of the phase
        sdphase = sqrt( (2*jksamples-2)*(1-abs(mean(phasefactorA,4))) );
    end
    
    if options.rlowess,
        for j=1:ndim,
            psd(:,j,j) = mslowess(f', psd(:,j,j));
            psderr(1,:,j,j) = mslowess(f', squeeze(psderr(1,:,j,j))');
            psderr(2,:,j,j) = mslowess(f', squeeze(psderr(2,:,j,j))');
            for l=1:ndim,
                coh(:,j,l) = mslowess(f', coh(:,j,l));
                coherr(1,:,j,l) = mslowess(f', squeeze(coherr(1,:,j,l))');
                coherr(2,:,j,l) = mslowess(f', squeeze(coherr(2,:,j,l))');
                pdc(:,j,l) = mslowess(f', pdc(:,j,l));
                pdcerr(1,:,j,l) = mslowess(f', squeeze(pdcerr(1,:,j,l))');
                pdcerr(2,:,j,l) = mslowess(f', squeeze(pdcerr(2,:,j,l))');
            end
        end
    end
    
    fit.state(k).psd = psd;
    fit.state(k).psderr = psderr;
    fit.state(k).coh = coh;
    fit.state(k).coherr = coherr;
    fit.state(k).pdc = pdc;
    fit.state(k).pdcerr = pdcerr;    
    fit.state(k).phase = phase;
    fit.state(k).sdphase = sdphase;  
    fit.state(k).f = f;  
end

end

%-------------------------------------------------------------------

function [tapers,eigs]=dpsschk(tapers,N,Fs)
% calculates tapers and, if precalculated tapers are supplied,
% checks that they (the precalculated tapers) the same length in time as
% the time series being studied. The length of the time series is specified
% as the second input argument N. Thus if precalculated tapers have
% dimensions [N1 K], we require that N1=N.
%
% From Chronux

if nargin < 3; error('Need all arguments'); end
sz=size(tapers);
if sz(1)==1 && sz(2)==2;
    [tapers,eigs]=dpss(N,tapers(1),tapers(2));
    tapers = tapers*sqrt(Fs);
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;
end

%-------------------------------------------------------------------

function [f,findx]=getfgrid(Fs,nfft,fpass)
% gets the frequency grid associated with a given fft based computation
% Called by spectral estimation routines to generate the frequency axes
%
% From Chronux

if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
    findx=find(f>=fpass(1) & f<=fpass(end));
else
    [~,findx]=min(abs(f-fpass));
end;
f=f(findx);
end

%-------------------------------------------------------------------

function J=mtfftc(data,tapers,nfft,Fs)
% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input:
%       data (in form samples x channels/trials or a single vector)
%       tapers (precalculated tapers from dpss)
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
%
% From Chronux

if nargin < 4; error('Need all input arguments'); end;
dtmp=[];
if isstruct(data);
    C=length(data);
    if C==1;
        fnames=fieldnames(data);
        eval(['dtmp=data.' fnames{1} ';'])
        data=dtmp(:);
    end
else
    [N,C]=size(data);
    if N==1 || C==1;
        data=data(:);
    end;
end;
[NC,C]=size(data); % size of data
[NK, K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers
J=fft(data_proj,nfft)/Fs;   % fft of projected data
end

%-------------------------------------------------------------------
function pdc = subrutpdc(S,numIterations,tol)
% obtains approximated pdc and dc from cross-spectral information

Nf = size(S,1); N = size(S,2);
% Wilson factorization
H = wilsonfact(permute(S,[2 3 1]),numIterations,tol);
% Computing the PDC
G = zeros(size(H));
pdc = zeros(size(H));
for f=1:Nf
    G(:,:,f) = inv(H(:,:,f));
    for i=1:N
        for j=1:N
            if i~=j
                pdc(i,j,f) = abs(G(i,j,f)) / sqrt(sum( abs(G(:,j,f)).^2 )) ;
            end
        end
    end
end
pdc = permute(pdc,[3 1 2]);
end

%-------------------------------------------------------------------
function [H, Z, S] = wilsonfact(S,Niterations,tol)
%
% This function is an implemention of Wilson's algorithm (Eq. 3.1)
% for spectral matrix factorization
%
% Inputs : S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency)
%        : fs (sampling frequency in Hz)
%        : freq (a vector of frequencies) at which S is given
% Outputs: H (transfer function)
%        : Z (noise covariance)
%        : psi (left spectral factor)
%
% Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities,"
% SIAM J. Appl. Math.23,420-426(1972).
% Modification over the function sfactorization_wilson (fieldtrip),
% implemented by  M. Dhamala & G. Rangarajan, UF, Aug 3-4, 2006.

% number of channels
m   = size(S,1);
N   = size(S,3)-1;
N2  = 2*N;

% preallocate memory for efficiency
Sarr   = zeros(m,m,N2) + 1i.*zeros(m,m,N2);
gam    = zeros(m,m,N2);
gamtmp = zeros(m,m,N2);
psi    = zeros(m,m,N2);
I      = eye(m); % Defining m x m identity matrix

%Step 1: Forming 2-sided spectral densities for ifft routine in matlab
for f_ind = 1:N+1
    Sarr(:,:,f_ind) = S(:,:,f_ind);
    if(f_ind>1)
        Sarr(:,:,2*N+2-f_ind) = S(:,:,f_ind).';
    end
end

%Step 2: Computing covariance matrices
for k1 = 1:m
    for k2 = 1:m
        gam(k1,k2,:) = real(ifft(squeeze(Sarr(k1,k2,:))));
    end
end

%Step 3: Initializing for iterations
gam0 = gam(:,:,1);
[h, dum] = chol(gam0);
if dum
    warning('initialization for iterations did not work well, using arbitrary starting condition');
    h = rand(m,m); h = triu(h); %arbitrary initial condition
end

for ind = 1:N2
    psi(:,:,ind) = h;
end

%Step 4: Iterating to get spectral factors
for iter = 1:Niterations
    for ind = 1:N2
        invpsi     = inv(psi(:,:,ind));
        g(:,:,ind) = invpsi*Sarr(:,:,ind)*invpsi'+I;%Eq 3.1
    end
    gp = PlusOperator(g,m,N+1); %gp constitutes positive and half of zero lags
    psi_old = psi;
    leaveloop = 0;
    for k = 1:N2
        psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
        %if isnan(rcond(psi(:,:,k))),
        %   leaveloop=1;
        %   break
        %end
        psierr(k)  = norm(psi(:,:,k)-psi_old(:,:,k),1);
    end
    %if leaveloop,
    %    psi = psi_old;
    %    break
    %end
    psierrf = mean(psierr);
    if(psierrf<tol),
        break;
    end; % checking convergence
end

%Step 5: Getting covariance matrix from spectral factors
for k1 = 1:m
    for k2 = 1:m
        gamtmp(k1,k2,:) = real(ifft(squeeze(psi(k1,k2,:))));
    end
end

%Step 6: Getting noise covariance & transfer function (see Example pp. 424)
A0    = gamtmp(:,:,1);
A0inv = inv(A0);
Z     = A0*A0.'; %Noise covariance matrix not multiplied by sampling frequency

H = zeros(m,m,N+1) + 1i*zeros(m,m,N+1);
for k = 1:N+1
    H(:,:,k) = psi(:,:,k)*A0inv;       %Transfer function
    %S(:,:,k) = psi(:,:,k)*psi(:,:,k)'; %Updated cross-spectral density
end
end

%---------------------------------------------------------------------
function gp = PlusOperator(g,nchan,nfreq)

g   = transpose(reshape(g, [nchan^2 2*(nfreq-1)]));
gam = ifft(g);

% taking only the positive lags and half of the zero lag
gamp  = gam;
beta0 = 0.5*gam(1,:);

gamp(1,          :) = reshape(triu(reshape(beta0, [nchan nchan])),[1 nchan^2]);
gamp(nfreq+1:end,:) = 0;

% reconstituting
gp = fft(gamp);
gp = reshape(transpose(gp), [nchan nchan 2*(nfreq-1)]);
end
