function [C,M] = osl_cov(D)
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
% [C,M] = osl_cov(D) also returns the mean
%
% Adam Baker 2014


%{
% PROTOTYPE CODE
% FOR FUTURE WORK
%
% The idea is that you can compute source space voxel covariance by getting 
% the sensor covariance, C_s, and then transforming it by tra. That is, for M_s and C_s in source space,
% the sensor space corresponding outputs M and C are
% M = tra*M_s
% C = tra*C_s*tra'
% This is way faster and will do away with memory issues in source space. But implementation should ideally be done when
% this function is encountered in production

function [C,M] = osl_cov(D,chaninds)
   if nargin < 2 || isempty(chaninds) 
       chaninds = [];
   end
   
   samples2use = good_samples(D,chaninds); % Check all good channels by default - this will be sensible in source space. In sensor space, probably better to just compute cov directly?


   if D.montage('getindex')
        montaged=true
       mont = D.montage('getmontage');
       D = D.montage('switch',0); 
       M = mont.tra*q2;
       else
        montaged = false
   end

   dat = D(:,:,:); % Read in sensor data
   M = zeros(D.nchannels,D.ntrials,D.nfrequencies);
   C = zeros(D.nchannels,D.nchannels,D.ntrials,D.nfrequencies);



   for j = 1:nfreqs
       for k = 1:ntrials
            if any(samples2use(:,k))
                if isa(D,'meeg') && isequal(D.transformtype,'TF')
                    C(:,:,k,f) = cov(dat(:,:,k,j)) % Or something like this, depends on what a TF MEEG looks like??
                    M() = mean()
                else
                    C(:,:,k,f) = cov(dat(:,:,k)) 
                    M() = mean()
                end
            end
        end
    end

    if montaged
        % Apply tra
        % Need to loop this over frequencies and trials
        M = mont.tra*M
        C=mont.tra*C*mont.tra';
    end
%}


if isa(D,'meeg')
    nfreqs   = max([1,D.nfrequencies]);
    nchans   = D.nchannels;
    ntrials  = D.ntrials;
    samples2use = good_samples(D);
else
    [nchans,nsamples] = size(D);
    ntrials = 1;
    nfreqs  = 1;
    if nchans > nsamples
        warning(['Input has ' num2str(nsamples) ' rows and ' num2str(nchans) ' columns. Consider transposing']);
    end
    samples2use = true(1,nsamples,1);
end

M = zeros(nchans,ntrials,nfreqs);
C = zeros(nchans,nchans,ntrials,nfreqs);
   
for f = 1:nfreqs
             
    for trl = 1:ntrials
        nsamples = sum(samples2use(1,:,trl));
        
        if nsamples,            
            % Compute means
            chan_blks = osl_memblocks([nchans,sum(samples2use(1,:,trl))],1);
            for i = 1:size(chan_blks,1)
                if isa(D,'meeg') && isequal(D.transformtype,'TF')
                    Dblk = D(chan_blks(i,1):chan_blks(i,2),f,:,trl);
                    Dblk = Dblk(:,:,samples2use(1,:,trl));
                else
                    Dblk = D(chan_blks(i,1):chan_blks(i,2),:,trl);
                    Dblk = Dblk(:,samples2use(1,:,trl));
                end
                Dblk = squeeze(Dblk);
                M(chan_blks(i,1):chan_blks(i,2),trl,f) = mean(Dblk,2);
            end

            % Compute covariance
            smpl_blks = osl_memblocks([nchans,nsamples],2);
            for i = 1:size(smpl_blks,1)                
                samples2use_blk = find(samples2use(1,:,trl));
                samples2use_blk = samples2use_blk(1,smpl_blks(i,1):smpl_blks(i,2));

                if isa(D,'meeg') && isequal(D.transformtype,'TF')
                    Dblk = squeeze(D(:,f,:,trl));
                    Dblk = Dblk(:,samples2use_blk);
                else
                    Dblk = D(:,:,trl);
                    Dblk = Dblk(:,samples2use_blk);
                end

                Dblk = bsxfun(@minus,Dblk,M(:,trl,f));
                C(:,:,trl,f) = C(:,:,trl,f) + Dblk*Dblk';
            end

            C(:,:,trl,f) = C(:,:,trl,f)./(nsamples-1);
        else % trial is bad or empty
            C(:,:,trl,f) = zeros(nchans);
        end
        
    end
end

end










