function W = transform_dipole_orientation(M, V, vectorBase, vectorLength)
%TRANSFORM_DIPOLE_ORIENTATION applies affine transformation to vector
%
% W = TRANSFORM_DIPOLE_ORIENTATION(M, V, VECTORBASE, VECTORLENGTH)
%   transforms vectors in the rows of V to W by applying affine
%   transformation matrix M. 
%
%   The transformation is effected by transforming the start and end points
%   of the vector, with start point set by VECTORBASE. The length of the
%   vector is given in VECTORLENGTH. If this is less than the granularity
%   of the simulation, or the grid spacing, an accurate transformation is
%   expected. 
% 
% See also: INVERSE_TRANSFORM_DIPOLE_ORIENTATION. 


%   Copyright 2014 OHBA
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


%   $LastChangedBy: giles.colclough@gmail.com $
%   $Revision: 213 $
%   $LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%   Contact: giles.colclough 'at' eng.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 25-Nov-2013 11:49:41

%% Input parsing
nVecs = size(V,1);

assert(size(V,2) == 3, ...
       [mfilename ':UnrecognisedVectorInput'], ...
       'Please enter 3-component vectors in rows of V. \n');
   
if size(vectorBase,1) == 1,
    vectorBase = repmat(vectorBase, nVecs, 1);
elseif size(vectorBase,1) ~= nVecs,
    error([mfilename 'VectorAndBaseMismatch'], ...
          'Please provide as many vector bases as vectors. \n');
end%if

if numel(vectorLength) == 1,
    vectorLength = repmat(vectorLength, nVecs, 1);
elseif numel(vectorBase) ~= nVecs,
    error([mfilename 'VectorAndLengthMismatch'], ...
          'Please provide as many vector lengths as vectors. \n');
end%if

%% Transform vector
% ensure vectors are unit length
VV = normalise_vectors(V, 2);

% define start and end points
vectorTip = vectorBase + scale_rows(VV, vectorLength);

% transform base and tip
newBase = spm_eeg_inv_transform_points(M, vectorBase);
newTip  = spm_eeg_inv_transform_points(M, vectorTip);

% reconstruct vectors
newV = newTip - newBase;

% normalise
W = normalise_vectors(newV, 2);
end%transform_dipole_orientation

%%
function y = scale_rows(x,s)
% SCALE_ROWS      Scale each row of a matrix.
% SCALE_ROWS(x,s) returns matrix y, same size as x, such that
% y(i,:) = s(i)*x(i,:)
% It is more efficient than diag(s)*x.

y = repmat(s(:), 1, size(x,2)).*x;
end%scale_rows

%%
function VV = normalise_vectors(V, dim)
%NORMALISE_VECTORS normalises rows or columns of a matrix
%
% NORMV = NORMALISE_VECTORS(V, DIM) normalises vectors along dimension 
%   DIM of V. If DIM=1, this treats columns as vectors and if DIM=2, this
%   treats rows as vectors. 
%
% NORMV = NORMALISE_VECTORS(V) is the same as NORMALISE_VECTORS(V, 1)
%
% Example:
%   If X = [3 3 3]
%   Then NORMALISE_VECTORS(X,1) is [1 1 1] and NORMALISE_VECTORS(X,2)
%   is [0.5774 0.5774 0.5774]. 


%   Copyright 2014 OHBA
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


%   $LastChangedBy: giles.colclough@gmail.com $
%   $Revision: 213 $
%   $LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%   Contact: giles.colclough 'at' eng.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 25-Nov-2013 11:49:41

if nargin < 2 || ~exist('dim', 'var') || isempty(dim),
    dim = 1;
end%if

VV = bsxfun(@rdivide, V, sqrt(sum(V.^2, dim)));
end%normalise_vectors
% [EOF]