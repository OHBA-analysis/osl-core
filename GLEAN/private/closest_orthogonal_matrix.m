function [L, d, rho] = closest_orthogonal_matrix(A)
%CLOSEST_ORTHOGONAL_MATRIX  Computes closest orthogonal matrix
%
% L = CLOSEST_ORTHOGONAL_MATRIX(A) returns orthogonal matrix L which
%   is closest to A, as measured by the Frobenius norm of (L-A), such that
%   transpose(L) * L is non-negative diagonal. 
%
% [L, D, RHO] = CLOSEST_ORTHOGONAL_MATRIX(A) also returns the scaling
%   factors for each column of A, D, and the final converged square
%   distance between L and A, RHO. 

%	References:
%	R. Evans,
%	"http://empslocal.ex.ac.uk/people/staff/reverson/uploads/Site/procrustes.pdf"
%   "Orthogonal, but not Orthonormal, Procrustes Problems", 1997
%	
%   See also: symmetric_orthogonalise, svd. 

%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
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


%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 07-Aug-2014 13:27:52

% settings
MAX_ITER  = 2e2;
TOLERANCE = max(1, max(size(A)) * fast_svds(A,1)) * eps(class(A)); % the tolerance used for matrices within matlab
DEBUG     = false;
dInitial  = ones(cols(A), 1); % use closest orthonormal matrix to initialise. Alternative: dInitial = sqrt(diag(A' * A)); 

% tests for convergence
convergence = @(rho, prevRho) abs(rho - prevRho) <= TOLERANCE;

% declare memory and initialise
iter = 0;
d    = dInitial;
rho  = NaN(MAX_ITER, 1);
A_T  = ctranspose(A);

while iter < MAX_ITER,
    iter = iter + 1;
    
    % find orthonormal polar factor
    V = symmetric_orthogonalise(A * diag(d));
    
    % minimise rho = |A - V D|^2 w.r.t d
    d = diag(A_T * V);
    
    % new best orthogonal estimate
    L = V * diag(d);
    
    % compute error term
    E = A - L;
    rho(iter) = sum(diag( E' * E )); % inline trace
    
    if DEBUG,
        fprintf('%s: iter %4d \t rho %g \n', mfilename, iter, rho(iter));  %#ok<UNRCH>
    end
    
    if iter > 1 && convergence(rho(iter), rho(iter - 1)),
        break
    end%if
end%convergence loop

% tidy vector
rho(isnan(rho)) = [];

if isequal(iter, MAX_ITER),
    warning([mfilename ':MaxIterationsHit'],                              ...
            ['%s: hit max iterations: %d. ',                              ...
             'You may not have found the optimal orthogonal matrix. \n'], ...
             mfilename, MAX_ITER);
end%if

if DEBUG,
    figure('Name', 'Convergence path', 'Color', 'w');                      %#ok<UNRCH>
    semilogy(2:(length(rho)-1), rho(2:end-1) - rho(end), 'b+');
    xlabel('iteration');
    ylabel('\rho_i - \rho_{\inf}');
end%if

end%closest_orthogonal_matrix
% [EOF]

