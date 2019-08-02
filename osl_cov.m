function [C,M] = osl_cov( D, varargin )
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
% Adam Baker 2014, JH 2019

    if isa(D,'meeg')
        nfreqs = max([ 1, D.nfrequencies ]);
        [smask,cid,sid,tid] = good_samples(D,varargin{:});
        
        is_timefreq = strcmp(D.transformtype,'TF');
        if nfreqs > 1 && ~is_timefreq
            warning('There are several frequencies, but no time-frequency transform.');
        end
        
        nchans   = numel(cid);
        nsamples = numel(sid);
        ntrials  = numel(tid);
    else
        [nchans,nsamples] = size(D);
        ntrials = 1;
        nfreqs  = 1;
        
        is_timefreq = false;
        if nchans > nsamples
            warning(['Input has ' num2str(nsamples) ' rows and ' num2str(nchans) ' columns. Consider transposing']);
        end
        if nargin > 1
            warning('Additional arguments are ignored for matrix inputs.');
        end
        
        smask = true(1,nsamples,1);
        cid = 1:nchans;
        sid = 1:nsamples;
        tid = 1:ntrials;
    end
    
    M = zeros(nchans,ntrials,nfreqs);
    C = zeros(nchans,nchans,ntrials,nfreqs);

    for k = 1:ntrials
        
        t = tid(k);
        s = sid(smask(1,:,k));
        
        nsmpl = numel(s);
        if nsmpl == 0, continue; end
        
        chan_blks = osl_memblocks([nchans,nsmpl],1);
        smpl_blks = osl_memblocks([nchans,nsmpl],2);
        
        for f = 1:nfreqs
            
            % Compute means
            for i = 1:size(chan_blks,1)
                b = chan_blks(i,1) : chan_blks(i,2);
                c = cid(b);
                if is_timefreq
                    Dblk = squeeze(D(c,f,s,t));
                else
                    Dblk = D(c,s,t);
                end
                M(b,k,f) = mean(Dblk,2);
            end

            % Compute covariance
            if nsmpl <= 1, continue; end
            for i = 1:size(smpl_blks,1)
                b = smpl_blks(i,1) : smpl_blks(i,2);

                if is_timefreq
                    Dblk = squeeze(D(cid,f,s(b),t));
                else
                    Dblk = D(cid,s(b),t);
                end

                Dblk = bsxfun( @minus, Dblk, M(:,k,f) );
                C(:,:,k,f) = C(:,:,k,f) + Dblk*Dblk';
            end
            C(:,:,k,f) = C(:,:,k,f)./(nsmpl-1);
            
        end

    end

end

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
