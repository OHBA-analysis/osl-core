function dipoleOri = get_orientations_from_mesh(dipolePos, brainMesh)
%GET_ORIENTATIONS_FROM_MESH
% finds the closest point on the cortical surface to the specified points.
% Uses the mesh normal as an orientation for the dipole. 


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
%   Contact: giles.colclough 'at' magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45


% find closest point on cortical surface
dipMeshInd = MEGsim.find_nearest_coordinate(dipolePos, brainMesh.vert);

% find the mesh normals
meshNorms  = spm_mesh_normals(struct('faces', brainMesh.face, ...
                                     'vertices', brainMesh.vert), ...
                              true);
                            
dipoleOri  = meshNorms(dipMeshInd, :);
end%get_orientations_from_mesh