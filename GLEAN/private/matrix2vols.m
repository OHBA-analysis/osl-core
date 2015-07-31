function [dat4d]=matrix2vols(dat2d,mask3d)
%takes 2d matrix (space x time) and 3d mask (xyz)
%and returns 4d data (xyz x time)

mask2d = reshape(mask3d,numel(mask3d),1) > 0;
dat4d  = zeros(numel(mask3d),size(dat2d,2));
dat4d(mask2d,:) = dat2d;
dat4d = reshape(dat4d,size(mask3d,1),size(mask3d,2),size(mask3d,3),size(dat2d,2));