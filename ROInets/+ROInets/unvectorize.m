function M = unvectorize(v)
%UNVECTORIZE convert edges into square matrix
%
% M = UNVECTORIZE(E) converts edges E into square symmetric network matrix
%   M. Diagonal entries will be zero. 
 
if size(v,2)>1,
    for m = size(v,2):-1:1,
        M(:,:,m) = ROInets.unvectorize(v(:,m));
    end%for
    return
end%if

nNodes   = ROInets.num_nodes(length(v));
M        = zeros(nNodes);
triUpInd = triu(true(nNodes),1);

M(triUpInd) = v(:);
M           = M + M';
end%unvectorize