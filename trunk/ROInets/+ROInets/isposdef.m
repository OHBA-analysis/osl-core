function b = isposdef(a)
% ISPOSDEF   Test for positive definite matrix.
%    ISPOSDEF(A) returns 1 if A is positive definite, 0 otherwise.
%    Using chol is much more efficient than computing eigenvectors.

% Written by Tom Minka

[~,p] = chol(a);
b = (p == 0);
