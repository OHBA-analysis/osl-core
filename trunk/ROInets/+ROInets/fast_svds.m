function [U, S, V] = fast_svds(X, N)
%FAST_SVDS RAM/time-efficient version of SVDS, singular value decomposition
% 
% S = FAST_SVDS(X, N) computes N components of the singular value
%   decomposition of X: X = U * S * V' where S is diagonal, and returns the
%   singular values in S.
%
% [U, S, V] = FAST_SVDS(X, N) computes N components of the singular value
%   decomposition of X: X = U * S * V' where S is diagonal. 
%   If N <= 0, rank(X) - abs(N) components are estimated. 
%
%   Note: no demeaning of X takes place within this function

%	Copyright 2013-4 OHBA, FMRIB
%	This program is free software: you can redirstribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 230 $
%	$LastChangedDate: 2014-08-07 20:52:19 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written by Steve Smith, FMRIB, 2013.

% parse component selection
if N < 1
  N = max(min(size(X)) + N, 1);
end%if

% Compute svd from eigenvalue decomposition of X'*X. 
% Choose which decomposition to take to maximise efficiency. 
if ROInets.rows(X) < ROInets.cols(X),
  if N < ROInets.rows(X),
      [U, d] = eigs(X * X', N);
      
  else
      [U, d] = eig(X * X'); 
      U      = fliplr(U); 
      d      = rot90(d, 2);
  end%if
  
  S = sqrt(abs(d));
  V = X' * (U * diag((1.0 / diag(S)))); 

else
  if N < ROInets.cols(X),
    [V, d] = eigs(X' * X, N);
    
  else
    [V, d] = eig(X' * X); 
    V      = fliplr(V); 
    d      = rot90(d, 2);
  end%if
  
  S = sqrt(abs(d));
  U = X * (V * diag((1.0 / diag(S)))); 

end%if

% change output based on required behaviour
if 1 == nargout, 
    U = S;
end%if
end%fast_svds
% [EOF]