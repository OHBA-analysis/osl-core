function n = num_nodes(nEdges)
%NUM_NODES returns the number of nodes in a network based on the number of
%edges
%
% N = NUM_NODES(NEDGES) returns the number of nodes (variables) in a square
%   symmetric network matrix with NEDGES upper off-diagonal elements
n = (sqrt(8*nEdges + 1) + 1) ./ 2;
end%num_nodes