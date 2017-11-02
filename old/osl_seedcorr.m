function C = osl_seedcorr(D,seedind,targetind,samples2use)
% Computes seed based correlation between a seed and a number of targets. 
% Note that this function does not perfrom any enveloping functionality
% (see e.g. osl_hilbenv.m)
%
% D             - MEEG object containing the data
% seedind       - index into channel/voxel in D for the seed
% targetind     - indices into channel(s)/voxel(s) for the target
% samples2use   - optionally specify samples to use (logical array)
%
% Adam Baker 2014

% Check 
try prcorr(1,1);
    useprcorr = 1;
catch
    useprcorr = 0;
end

if nargin < 4
    samples2use = true(1,D.nsamples);
end

samples2use = samples2use(:)';

samples2use = samples2use & goodsamples(D);
samples2use = find(samples2use);

C = zeros(length(targetind),D.ntrials);
for trl = 1:D.ntrials
    
    % Read data from file in blocks
    blks = osl_memblocks(size(D),1);
    
    seed_data = squeeze(D(seedind,samples2use,trl))';
    
    for iblk = 1:size(blks,1)
        
        dat_blk = D(targetind(blks(iblk,1):blks(iblk,2)),samples2use,trl)';
        
        for i = 1:size(dat_blk,2)
            if useprcorr
                C(blks(iblk,1)+i-1,trl) = prcorr(seed_data,dat_blk(:,i));
            else
                C(blks(iblk,1)+i-1,trl) = corr(seed_data,dat_blk(:,i));
            end
        end
    end
end


end
