function [ index, vec, dist ] = nearest_vec( vec_array, vec_to_find )

% [ index, vec, dist ] = nearest_vec( vec_array, vec_to_find )
%
% MWW

dists=(sqrt(sum((vec_array-repmat(vec_to_find,size(vec_array,1),1)).^2,2)));

[dist,index]=min(dists);

vec=vec_array(index,:);

end

