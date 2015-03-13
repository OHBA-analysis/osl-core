function C = osl_cov(X)
% Computes the covariance of a [channels x samples] matrix without 
% encountering memory issues. This allows the covariance matrix to be 
% computed for large data matrices without running out of memory.
% 
% This function can also compute the (channels x channels) covariance 
% matrix of a (channelx x samples x trials) MEEG object (using only good 
% samples). 
%
% Usage:
% C = osl_cov(X)
% 
% OR:
% C = osl_cov(D)
%
% Adam Baker 2014

if isa(X,'meeg')
    nch  = X.nchannels;
    %samples2use = reshape(squeeze(~all(badsamples(X,':',':',':'))),X.nsamples*X.ntrials,1); 
    samples2use = squeeze(~all(badsamples(X,':',':',':')));
    samples2use = samples2use(:);
    nsmp = sum(samples2use);
else
    [nch,nsmp] = size(X);
    if nch > nsmp
        error(['Input has ' num2str(nsmp) ' rows and ' num2str(ncol) ' columns. Consider transposing']);
    end
    samples2use = true(nsmp,1);
end

% Compute means
chan_blks = osl_memblocks(X,1);
M = zeros(nch,1);
for i = 1:size(chan_blks,1)
    Xblk = X(chan_blks(i,1):chan_blks(i,2),:,:);
    Xblk = Xblk(:,samples2use);
    M(chan_blks(i,1):chan_blks(i,2)) = mean(Xblk,2);
end

% Compute covariance
smpl_blks = osl_memblocks(X,2);
C = zeros(nch);
for i = 1:size(smpl_blks,1)
    Xblk = X(:,smpl_blks(i,1):smpl_blks(i,2),:);
    Xblk = bsxfun(@minus,Xblk,M);
    samples2use_blk = samples2use(smpl_blks(i,1):smpl_blks(i,2),:);
    C = C + Xblk(:,samples2use_blk)*Xblk(:,samples2use_blk)';
end
C = C./(nsmp-1);


end