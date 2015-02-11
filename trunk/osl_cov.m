function C = osl_cov(X, divFlag)
% Computes the covariance of a [channels x samples] matrix without 
% encountering memory issues. This allows the covariance matrix to be 
% computed for large data matrices without running out of memory.
% 
% This function can also compute the (channels x channels) covariance 
% matrix of a (channelx x samples x trials) MEEG object (using only good 
% samples). 
%
% Usage:
% C = osl_cov(X) computes the covariance of matrix X. 
%
% C = osl_cov(D) computes the covariance of good MEG channels in D. 
%
% C = osl_cov(..., 1) computes the covariance using biased estimator
%   X*X'/n, instead of dividing by (n-1). 
%
% Adam Baker 2014

% parse choice of divisor
if nargin < 2 || ~exist('divFlag', 'var')
    divFlag = 0;
end%if

% parse behaviour based on input
if isa(X,'meeg')
    samples2use = reshape(squeeze(~all(badsamples(X,':',':',':'))),X.nsamples*X.ntrials,1);
    nsmp        = sum(samples2use);
    megchannels = X.indchantype({'MEG', 'MEGANY'}, 'GOOD'); % need both MEG and MEGANY to cover all bases.
    nch         = length(megchannels);
else
    [nch,nsmp] = size(X);
    if nch > nsmp
        error(['Input has ' num2str(nsmp) ' rows and ' num2str(ncol) ' columns. Consider transposing']);
    end
    samples2use = true(1,nsmp);
    megchannels = 1:nch;
end

% Compute means
chan_blks = osl_memblocks(X,1);
M = zeros(nch,1);
for i = 1:size(chan_blks,1)
    % ensure we only look at good meg channels
    blk_chan_inds = intersect(chan_blks(i,1):chan_blks(i,2), megchannels);
    
    Xblk = X(blk_chan_inds,:,:);
    Xblk = reshape(Xblk,size(Xblk,1),[]);
    
    M(blk_chan_inds) = mean(Xblk(:,samples2use), 2);
end

% Compute covariance
smpl_blks = osl_memblocks(X, 2);
C         = zeros(nch);
for i = 1:size(smpl_blks,1)
    Xblk = X(megchannels, smpl_blks(i,1):smpl_blks(i,2), :);
    Xblk = reshape(Xblk, size(Xblk,1), []);
    Xblk = bsxfun(@minus, Xblk, M);
    
    samples2use_blk = samples2use(smpl_blks(i,1):smpl_blks(i,2));
    
    C = C + Xblk(:,samples2use_blk) * Xblk(:,samples2use_blk)';
end

if divFlag, % biased estimate
    C = C ./ nsmp;
else % standard unbiased estimator
    C = C ./ (nsmp-1);
end%if
end