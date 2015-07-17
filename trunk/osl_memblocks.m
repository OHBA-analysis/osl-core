function blks = osl_memblocks(SZE,dim,MaxMem)
% Returns indices into matrix of size SZE that splits the data into blocks no greater than 
% MAXMEM bytes along dimension DIM

if nargin < 2
    dim = 1;
end

if nargin < 3
    MaxMem = 200*2^20; % bytes
end

Mem_blk = 8*prod(SZE(setxor(dim,1:length(SZE)))); %#ok

blk_size = floor(MaxMem./Mem_blk);
blks     = 0:blk_size:SZE(dim);
blks     = unique([blks SZE(dim)]);
blks     = [blks(1:end-1)+1; blks(2:end)]';

end