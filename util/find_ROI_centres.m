function coords = find_ROI_centres(spatialMap, spatialRes, isBinary, OSLDIR)
%FIND_ROI_CENTRES  finds centre of mass of a set of ROIs


%   Copyright 2014 Giles Colclough
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


%   $LastChangedBy$
%   $Revision$
%   $LastChangedDate$
%   Contact: giles.colclough@eng.ox.ac.uk
%   Originally written on: MACI64 by Giles Colclough, 29-Jul-2014 11:39:10

if nargin<4
    global OSLDIR 
end

[nVoxels, nParcels] = size(spatialMap);

brainMaskName = fullfile(OSLDIR, 'std_masks', sprintf('MNI152_T1_%dmm_brain_mask.nii.gz', spatialRes));
MNIcoords     = osl_mnimask2mnicoords(brainMaskName);
assert(ROInets.rows(MNIcoords) == nVoxels);

for iParcel = nParcels:-1:1, 
    map = spatialMap(:, iParcel);
    
    % find ROI
    if isBinary,
        cutOff = 0;        
    else
        % extract top 5% of values
        cutOff = prctile(map, 95);
    end%if
    ROIinds = (map > cutOff);
    
    % find weightings
    if isBinary,
        masses = ones(sum(ROIinds), 1);
    else
        masses = map(ROIinds);
    end
    
    % find CoM
    CoM = (masses' * MNIcoords(ROIinds,:)) ./ sum(masses);
    
    coords(iParcel, :) = CoM;
end%for
end%find_ROI_centres