function L = symmetric_orthogonalise(A, maintainMagnitudes)
%SYMMETRIC_ORTHOGONALISE closest orthogonal matrix
% 
% L = SYMMETRIC_ORTHOGONALISE(A) returns orthonormal matrix L which
%   is closest to A, as measured by the Frobenius norm of (L-A). 
%
%   The orthogonal matrix is constructed from a singular value decomposition
%   of A. 
%
% L = SYMMETRIC_ORTHOGONALISE(A, KEEP_MAGNITUDES) returns the orthogonal
%   matrix L, whose columns have the same magnitude as the respective
%   columns of A, and which is closest to A, as measured by the Frobenius
%   norm of (L-A), if KEEP_MAGNITUDES is TRUE. 
%
%   The orthogonal matrix is constructed from a singular value decomposition
%   of A. 
%
%   See also: ROINETS.HOUSEHOLDER_ORTHOGONALISE, ORTH, SVD. 

% References: Naidu, A. R. "Centrality of Lowdin Orthogonalizations",
%   arXiv 1105.3571v1, May 2011. 
%   Available at: http://arxiv-web3.library.cornell.edu/pdf/1105.3571v1.pdf

%	Copyright 2013 OHBA
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
%	$Revision: 239 $
%	$LastChangedDate: 2014-08-15 14:58:49 +0100 (Fri, 15 Aug 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 31-Oct-2013 13:30:05

if nargin < 2 || ~exist('maintainMagnitudes', 'var'),
    maintainMagnitudes = false;
end

if maintainMagnitudes,
    D = diag(sqrt(diag(A' * A)));
    
    % call function again
    Lnorm = ROInets.symmetric_orthogonalise(A * D, false);
    
    % scale result
    L = Lnorm * D;
    
else
    [U, S, V] = svd(A, 'econ');

    if ~isempty(S),
        % we need to check that we have sufficient rank
        S   = diag(S);
        tol = max(size(A)) * S(1) * eps(class(A));
        r   = sum(S > tol);
        
        isFullRank = (r >= ROInets.cols(A));

        if isFullRank,
            % polar factors of A
            L = U * conj(transpose(V));

        else % not enough degrees of freedom
            error([mfilename ':RankError'], ...
                  ['The input matrix is not full rank. \n', ...
                   '    There are not enough degrees of freedom to perform a ', ...
                   'sensible orthogonalisation. \n', ...
                   '    Try a different orthogonalisation method. \n']);
        end%if        

    else
        L = [];
    end%if
end%if
end%symmetric_orthogonalise
% [EOF]