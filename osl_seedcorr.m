function C = osl_seedcorr(D,seedind,targetind)
% Computes seed based correlation between a seed and a number of targets. 
% Note that this function does not perfrom any enveloping functionality
% (see e.g. osl_hilbenv.m)
%
% D         - MEEG object containing the data
% seedind   - index into channel/voxel in D for the seed
% targetind - indices into channel(s)/voxel(s) for the target
%
% Adam Baker 2014

% Check 
try prcorr(1,1);
    useprcorr = 1;
catch
    useprcorr = 0;
end

C = zeros(length(targetind),D.ntrials);
for trl = 1:D.ntrials
    
    % Read data from file in blocks
    Mem_max = 200*2^20;
    Mem_chan = 8*numel(D(1,:,trl));
    blk_size = floor(Mem_max./Mem_chan);
    blks = 1:blk_size:length(targetind);
    blks = unique([blks length(targetind)]);
    blks = [blks(1:end-1); blks(2:end)]';
    
    seed_data = squeeze(D(seedind,:,trl))';
    
    for iblk = 1:size(blks,1)
        
        dat_blk = D(targetind(blks(iblk,1):blks(iblk,2)),:,trl)';
        
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
