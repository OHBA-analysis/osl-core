function [X, tol, r, s2] = mwpinv(A,varargin)
%MWPINV   Pseudoinverse with reduced rank
%   X = MWPINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS(class(A)).
%
%   MWPINV(A, R) calculates PINV(A) using a rank R for matrix A; i.e. using
%   only the first R singular values. 
% 
% 
%   Class support for input A: 
%      float: double, single
%
%   See also PINV, RANK.

% Written by Mark Woolrich, based on PINV by TMW. 
% mark.woolrich 'at' ohba.ox.ac.uk

if isempty(A)     % quick return
  X = zeros(size(A'),class(A));  
  return  
end

% [U,S,V] = svd(A,0); s = diag(S);
[m,n] = size(A);

if n > m
   X = pinv(A',varargin{:})';
else
   [U,S,V] = svd(A,0);
   if m > 1, s = diag(S);
      elseif m == 1, s = S(1);
      else s = 0;
   end
   if nargin == 2,
       
       r = varargin{1}; % rank to use
       tol=[];
       if r<0        

         diffs=abs(diff(s)./s(1:end-1));
         for i=2:length(diffs),
             if(diffs(i) > 6*mean(diffs(1:i))), 
                 break;
             end;
         end; 
         r=i;
       
       end;
      
   else
      %tol_old = max(m,n) * eps(max(s))
   
      tol=max(m,n) * eps(max(s));
   
      r = sum(s > tol);       % calculate rank using tolerance
      
   end
 
   s2=s;
   
   if (r == 0)
      X = zeros(size(A'),class(A));
   else
      s = diag(ones(r,1)./s(1:r));
      X = V(:,1:r)*s*U(:,1:r)';
   end
end
