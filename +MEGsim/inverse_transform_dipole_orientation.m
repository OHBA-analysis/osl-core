function W = inverse_transform_dipole_orientation(M, V, vectorBase, vectorLength)
%INVERSE_TRANSFORM_DIPOLE_ORIENTATION applies affine transformation in
% reverse to vector
%
% W = INVERSE_TRANSFORM_DIPOLE_ORIENTATION(M, V, VECTORBASE, VECTORLENGTH)
%   transforms vectors in the rows of V to W by applying affine
%   transformation matrix M in reverse. 
%
%   The transformation is effected by transforming the start and end points
%   of the vector, with start point set by VECTORBASE. The length of the
%   vector is given in VECTORLENGTH. If this is less than the granularity
%   of the simulation, or the grid spacing, an accurate transformation is
%   expected. 
% 
% See also: TRANSFORM_DIPOLE_ORIENTATION. 


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

Minv = MEGsim.invert_affine_transformation(M);

W = MEGsim.transform_dipole_orientation(Minv, V, vectorBase, vectorLength);
end%inverse_transform_dipole_orientation
% [EOF]
%{
%%% TEST SCRIPT %%%
% F = [-0.0377   -1.1565    0.0230    0.0368
%       1.0374   -0.0222   -0.2820  -13.2476
%       0.3861    0.0115    1.1625  -57.2605
%            0         0         0    1.0000];
%        
% vec     = [1 1 1]/sqrt(3)
% vecBase = [10 -12 22];
% L = 0.5;
%
% W = inverse_transform_dipole_orientation(F, vec, vecBase, L)
% newBase = inverse_affine_transform_points(F, vecBase);
% 
% shouldBeSameAsVec = transform_dipole_orientation(F, W, newBase, L)
%}