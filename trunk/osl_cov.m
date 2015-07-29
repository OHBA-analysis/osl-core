function C = osl_cov(D)
% Computes the covariance of a [channels x samples] matrix without
% encountering memory issues. This allows the covariance matrix to be
% computed for large data matrices without running out of memory.
%
% This function can also compute the (channels x channels x trials) covariance
% matrix of a (channelx x samples x trials) MEEG object (using only good
% samples).
%
% Usage:
% C = osl_cov(D)
%
% OR:
% C = osl_cov(D)
%
% Adam Baker 2014

if isa(D,'meeg')
    samples2use = permute(~all(badsamples(D,':',':',':')),[2 3 1]);
    nfreqs   = max([1,D.nfrequencies]);
    nchans   = D.nchannels;
    ntrials  = D.ntrials;
    if isrow(samples2use), samples2use = samples2use(:); end%if
else
    [nchans,nsamples] = size(D);
    ntrials = 1;
    nfreqs  = 1;
    if nchans > nsamples
        error(['Input has ' num2str(nsamples) ' rows and ' num2str(nchans) ' columns. Consider transposing']);
    end
    samples2use = true(nsamples,1);
end


M = zeros(nchans,ntrials,nfreqs);
C = zeros(nchans,nchans,ntrials,nfreqs);
   
for f = 1:nfreqs
             
    for trl = 1:ntrials
            
        nsamples = sum(samples2use(:,trl));
        
        % Compute means
        chan_blks = osl_memblocks([nchans,sum(samples2use(:,trl))],1);
        for i = 1:size(chan_blks,1)
            if isa(D,'meeg') && isequal(D.transformtype,'TF')
                Dblk = D(chan_blks(i,1):chan_blks(i,2),f,find(samples2use(:,trl)),trl);
            else
                Dblk = D(chan_blks(i,1):chan_blks(i,2),find(samples2use(:,trl)),trl);
            end
            Dblk = squeeze(Dblk);
            M(chan_blks(i,1):chan_blks(i,2),trl,f) = mean(Dblk,2);
        end
        
        % Compute covariance
        smpl_blks = osl_memblocks([nchans,nsamples],2);
        for i = 1:size(smpl_blks,1)
            
            samples2use_blk = find(samples2use(:,trl));
            samples2use_blk = samples2use_blk(smpl_blks(i,1):smpl_blks(i,2));
            
            if isa(D,'meeg') && isequal(D.transformtype,'TF')
                Dblk = squeeze(D(:,f,samples2use_blk,trl));
            else
                Dblk = D(:,samples2use_blk,trl);
            end
            
            Dblk = bsxfun(@minus,Dblk,M(:,trl,f));
            C(:,:,trl,f) = C(:,:,trl,f) + Dblk*Dblk';
        end
        C(:,:,trl,f) = C(:,:,trl,f)./(nsamples-1);
        
    end
end

end










