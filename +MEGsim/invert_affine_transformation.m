function Minv = invert_affine_transformation(M)
%INVERT_AFFINE_TRANSFORMATION inverts and affine transformation matrix
%
% MINV = INVERT_AFFINE_TRANSFORMATION(M) inverts affine transformation
%   matrix M. 


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

A = M(1:3,1:3);
b = M(1:3,4);

if det(A) == 0,
    error([mfilename ':InvalidInverse'], ...
          'Taking an inverse of this transform is impossible. \n');
end%if

Minv = [ pinv(A),    - pinv(A) * b; ...
         zeros(1,3),             1];
end%invert_affine_transformation
% [EOF]