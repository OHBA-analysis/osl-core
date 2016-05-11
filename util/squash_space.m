function [xr, mask] = squash_space(x, m)

% function [ret, mask] = squash_space(x)
%
% performs dims=size(x);ret = reshape(x,prod(dims(1:3)),dims(4));
% x needs to be 4D
%
% [ret, mask] = squash_space(x, m)
% masks and then squashes
% m is a 3D mask
%
% see squash
   
if(length(size(x))~=4)
    display('x needs to be a 4D matrix');
    return;
end;

[xr, dims]=sq_space(x);

if(nargin > 1)
   sqm = sq(m);
   coords = find(sqm>0);
   xr = xr(coords,:);
   mask = m;
else
    mask = ones(dims);
end;
   
function [xr, dims] = sq_space(x)
   
dims = size(x);
xr = reshape(x,prod(dims(1:3)),dims(4));

function [xr, dims] = sq(x)
   
dims = size(x);
xr = reshape(x,prod(dims),1);