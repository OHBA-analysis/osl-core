function fileNameOut = nii_parcel_quicksave(data, parcelFlag, parcelWeights, filename, varargin)
%NII_PARCEL_QUICKSAVE	Saves data in parcels as nifti
% NII_PARCEL_QUICKSAVE(DATA, PARCELFLAG, PARCELWEIGHTS, FILENAME, SPATIALRES) saves
%   DATA in nifti file FILENAME using spatial resolution SPATIALRES. The
%   DATA form an (nParcels) x (nVolumes) matrix and PARCELFLAG (nVoxels) x
%   (nParcels) is a binary matrix identifying the membership of voxels in
%   parcels, and PARCELWEIGHTS (nVoxels) x (nParcels) is a matrix
%   indicating the membership of each voxel in a parcel (set to an
%   empty matrix, to indicate a binary parcellation).
%
% NII_PARCEL_QUICKSAVE(DATA, PARCELFLAG, PARCELWEIGHTS, FILENAME, SPATIALRES, RESAMP, INTERP)
%   uses resampling method RESAMP and interpolation method INTERP in the
%   call to nii_quicksave. 
%
% See also: NII_QUICKSAVE. 


%	Copyright 2013 OHBA
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


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 213 $
%	$LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 10-Dec-2013 10:07:57

NII_MAX_SIZE = 32767;

% check save location exists
saveDir = fileparts(filename);
if ~exist(saveDir, 'dir'), 
    warning([mfilename ':CreatingSaveDirectory'], ...
            'Save location did not exist. Creating directory \n   %s\n', ...
            saveDir);
    ROInets.make_directory(saveDir); 
end%if

% strip extension from filename - actually, don't. This takes apart names
% with periods in the middle, but no explicit extension. 
% filename = fullfile(saveDir, fileStem);

% check parcelFlag is logical
assert(islogical(parcelFlag), ...
       [mfilename, ':nonLogicalInput'], ...
       'parcelFlag must be a logical array. \n');

% check that each voxel is only a member of one parcel
assert(~any(ROInets.row_sum(parcelFlag) > 1), ...
       [mfilename ':MultipleParcelOccupancy'], ...
       'Each voxel can be a member of at most one parcel. \n');
   
% check data sizes match
[nVoxels, nParcels] = size(parcelFlag);
assert(isequal(ROInets.rows(data), nParcels), ...
       [mfilename ':InputSizeMismatch'], ...
       'data must have rows equal to number of columns in parcelFlag. \n');

% declare memory to ensure correct size
rePackedData = zeros(nVoxels, ROInets.cols(data));

% repack data into voxel form
if isempty(parcelWeights)
    for iParcel = nParcels:-1:1,
        insertInds                  = parcelFlag(:, iParcel);
        rePackedData(insertInds, :) = repmat(data(iParcel, :), ...
                                             sum(insertInds), 1);
    end;
else
    rePackedData=parcelWeights*data;
end;

if ROInets.cols(rePackedData) < NII_MAX_SIZE,
    % save using osl function
    fileNameOut=nii_quicksave(rePackedData, filename, varargin{:});
    
else
    fprintf('Nii file limit exceeded, saving as .mat \n');
    fileNameOut = [filename '.mat']; % overwrite previous extension
    save(fileNameOut, 'rePackedData', '-v7.3');
end%if

end%nii_parcel_quicksave
% [EOF]