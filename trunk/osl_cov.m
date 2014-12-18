function C = osl_cov(X)
% Computes the covariance of a [channels x samples] matrix without 
% encountering memory issues. Can also compute the (channels x channels) 
% covariance matrix of a continuous MEEG object. This allows the covariance
% matrix to be computed for large data matrices without running out of
% memory.
%
% Usage:
% C = osl_cov(X)
%
%
% Adam Baker 2014

if isa(X,'meeg')
    nch = X.nchannels;
    nsmp = X.nsamples;
else
    [nch,nsmp] = size(X);
    if nch > nsmp
        error(['Input has ' num2str(nsmp) ' rows and ' num2str(ncol) ' columns. Consider transposing']);
    end
end

% Assess memory usage
Mem_max = 200*2^20;
Mem_smp = 8*nch;
Mem_ch = 8*nsmp;

% Compute means
blk_size = floor(Mem_max./Mem_ch);    
blks = 0:blk_size:nch;
blks = unique([blks nch]);
blks = [blks(1:end-1)+1; blks(2:end)]';
M = zeros(nch,1);
for i = 1:size(blks,1)
    Xblk = X(blks(i,1):blks(i,2),:);
    M(blks(i,1):blks(i,2)) = mean(Xblk,2);
end


% Compute covariance
blk_size = floor(Mem_max./Mem_smp);    
blks = 0:blk_size:nsmp;
blks = unique([blks nsmp]);
blks = [blks(1:end-1)+1; blks(2:end)]';

C = zeros(nch);
for i = 1:size(blks,1)
    Xblk = bsxfun(@minus,X(:,blks(i,1):blks(i,2)),M);
    C = C + Xblk*Xblk';
end
C = C./(nsmp-1);


end