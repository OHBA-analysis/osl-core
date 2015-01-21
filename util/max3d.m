function [a b c max_stat]=max3d(X)

% [indi indj indk]=max2d(X)

[tmp ind]=max(X(:));

[a b c]=ind2sub(size(X),ind);

max_stat=X(a,b,c);

end
