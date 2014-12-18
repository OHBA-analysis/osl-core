function fit = hmmspectramar(hmm,options)
% Get ML spectral estimates from MAR model
%
% INPUT
% X             time series 
% T             Number of time points for each time series
% hmm           An hmm-mar structure 
% options 

%  .loadings	In case data are principal components, the PCA loading matrix  
%  .Fs:       Sampling frequency
%  .fpass:    Frequency band to be used [fmin fmax] (default [0 fs/2])
%  .Nf        No. of frequencies to be computed in the range minHz-maxHz
%
% OUTPUT
% fit is a list with K elements, each of which contains: 
% fit.state(k).psd     [Nf x N x N] Power Spectral Density matrix
% fit.state(k).ipsd     [Nf x N x N] Inverse Power Spectral Density matrix
% fit.state(k).coh     [Nf x N x N] Coherence matrix
% fit.state(k).pcoh    [Nf x N x N] Partial Coherence matrix
% fit.state(k).pdc   [Nf x N x N] Baccala's Partial Directed Coherence
% fit.state(k).phase     [Nf x N x N] Phase matrix
% fit.state(k).f     [Nf x 1] Frequency vector
%
% Author: Diego Vidaurre, OHBA, University of Oxford (2014)

ndim = size(hmm.state(1).W.Mu_W,2); 
if ~isfield(options,'loadings'), options.loadings=eye(ndim); end;
if ~isfield(options,'Fs'), options.Fs=1; end;
if ~isfield(options,'fpass'),  options.fpass=[0 options.Fs/2]; end;
if ~isfield(options,'Nf'),  options.Nf=256; end;
loadings = options.loadings;
if hmm.train.whitening, loadings = loadings * iA; end

M = size(options.loadings,1);
freqs = (0:options.Nf-1)*( (options.fpass(2) - options.fpass(1)) / (options.Nf-1)) + options.fpass(1);
w=2*pi*freqs/options.Fs;

orders = 1:hmm.train.timelag:hmm.train.order; order = orders(end); 
K = length(hmm.state);
W = zeros(K,order,M,M);
for k=1:K
    for i=1:length(orders),
        W(k,i,:,:) = loadings *  hmm.state(k).W.Mu_W((1:ndim) + (i-1)*ndim,:) * loadings' ;
    end
end;

for k=1:K
    
    switch hmm.train.covtype
        case 'uniquediag'
            covmk = diag(hmm.Omega.Gam_rate / hmm.Omega.Gam_shape);
            preck=diag(hmm.Omega.Gam_shape ./ hmm.Omega.Gam_rate);
        case 'diag'
            covmk = diag(hmm.state(k).Omega.Gam_rate / hmm.state(k).Omega.Gam_shape);
            preck=diag(hmm.state(k).Omega.Gam_shape ./ hmm.state(k).Omega.Gam_rate);
        case 'uniquefull'
            covmk = hmm.Omega.Gam_rate ./ hmm.Omega.Gam_shape;
            preck= inv(covmk);
        case 'full'
            covmk = hmm.state(k).Omega.Gam_rate ./ hmm.state(k).Omega.Gam_shape;
            preck= inv(covmk);
    end
    
    % Get Power Spectral Density matrix and DTF for state K
    for ff=1:options.Nf,
        af_tmp=eye(ndim);
        for i=1:length(orders),
            o = orders(i);
            af_tmp=af_tmp - permute(W(k,i,:,:),[3 4 1 2]) * exp(-1i*w(ff)*o);
        end
        iaf_tmp=inv(af_tmp);
        fit.state(k).psd(ff,:,:) = iaf_tmp * covmk * iaf_tmp';
        fit.state(k).ipsd(ff,:,:) = inv(permute(fit.state(k).psd(ff,:,:),[3 2 1]));
        
        % Get PDC
        for j=1:ndim,
            prec_nn=1/sqrt(covmk(j,j));
            for l=1:ndim,
                fit.state(k).pdc(ff,j,l) = prec_nn * abs(af_tmp(j,l))/sqrt(abs(af_tmp(:,l)'*preck*af_tmp(:,l)));
            end
        end
    end
    
    % Get Coherence and Phase
    for j=1:ndim,
        for l=1:ndim,
            rkj=fit.state(k).psd(:,j,l)./(sqrt(fit.state(k).psd(:,j,j)).*sqrt(fit.state(k).psd(:,l,l)));
            %fit.state(k).coh(:,nn,n)=abs(rkj);
            fit.state(k).coh(:,j,l)=rkj;
            fit.state(k).pcoh(:,j,l)=-fit.state(k).ipsd(:,j,l)./(sqrt(fit.state(k).ipsd(:,j,j)).*sqrt(fit.state(k).ipsd(:,l,l)));
            fit.state(k).phase(:,j,l)=atan(imag(rkj)./real(rkj));
        end
    end
    
    for j=1:ndim,
        fit.state(k).psd(:,j,j) = real(fit.state(k).psd(:,j,j));
        fit.state(k).ipsd(:,j,j) = real(fit.state(k).ipsd(:,j,j));
    end
    
        fit.state(k).f = freqs;
end

