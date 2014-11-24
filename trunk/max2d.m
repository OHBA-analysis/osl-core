function [indi indj]=max2d(X)

% [indi indj]=max2d(X)

[v,ind]=max(X);
[v1,ind1]=max(max(X));
indi=ind(ind1);
indj=ind1;

end
