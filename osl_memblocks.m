function blks = osl_memblocks(X,dim,MaxMem)
% Returns indices into X that splits the data into blocks no greater than 
% MAXMEM bytes along dimension DIM

if nargin < 2
    dim = 1;
end

if nargin < 3
    MaxMem = 200*2^20; % bytes
end

Mem_blk = 8*prod(size(X,setxor(dim,1:length(size(X))))); %#ok

blk_size = floor(MaxMem./Mem_blk);
blks     = 0:blk_size:size(X,dim);
blks     = unique([blks size(X,dim)]);
blks     = [blks(1:end-1)+1; blks(2:end)]';

end