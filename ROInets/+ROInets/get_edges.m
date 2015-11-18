function e = get_edges(W)
%GET_EDGES retrieves edges from network matrix 
%
% E = GET_EDGES(W) retrieves the edges from square symmetric network matrix
%   W. 

% reshape W
[nEdges, ~, nSubs]  = size(W);
edgeInd = repmat(triu(true(nEdges),1), [1 1 nSubs]);
e       = reshape(W(edgeInd), [], nSubs);
end%get_edges