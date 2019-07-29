function [xr] = unsquash_space(x,m)

% function [x] = unsquash(x,m)
%
% [x] = unsquash(x,dims,m)
%
% unsquashes restoring background
% specified by m mask
%
% see squash

coords = find(squash(m)>0);
xr = zeros(prod(size(m)),size(x,2));
xr(coords,:) = x;

xr = reshape(xr,[size(m),size(x,2)]);

