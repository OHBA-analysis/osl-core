function fileNameOut = nii_parcel_quicksave(data, parcelFlag, filename, varargin)
%NII_PARCEL_QUICKSAVE	Saves data in parcels as nifti
% NII_PARCEL_QUICKSAVE(DATA, PARCELFLAG, FILENAME, OPTIONS) or
%   saves DATA in nifti file FILENAME using the spatial resolution specified in
%   OPTIONS.INPUT_SPAT_RES (in mm). See nii_quicksave for other OPTIONS. The
%   DATA form an (nParcels) x (nVolumes) matrix and PARCELFLAG (nVoxels) x
%   (nParcels) is a binary matrix identifying the membership of voxels in
%   parcels, or a matrix indicating the membership of each voxel in a
%   parcel.
% OR:
% NII_PARCEL_QUICKSAVE(DATA, PARCELFLAG, FILENAME, SPATIALRES, RESAMP, INTERP) 
%   Is the old interface. Where is the output spatial resolution specified
%   (in mm). See nii_quicksave for other settings
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
   
% check data sizes match
[nVoxels, nParcels] = size(parcelFlag);
assert(isequal(ROInets.rows(data), nParcels), ...
       [mfilename ':InputSizeMismatch'], ...
       'data must have rows equal to number of columns in parcelFlag. \n');

% declare memory to ensure correct size
rePackedData = zeros(nVoxels, ROInets.cols(data));

% repack data into voxel form
if islogical(parcelFlag),
    disp('Computing maps from binary parcellation');    
    
    % check that each voxel is only a member of one parcel
    assert(~any(ROInets.row_sum(parcelFlag) > 1), ...
       [mfilename ':MultipleParcelOccupancy'], ...
       'Each voxel can be a member of at most one parcel. \n');

    for iParcel = nParcels:-1:1,
        insertInds                  = parcelFlag(:, iParcel);
        rePackedData(insertInds, :) = repmat(data(iParcel, :), ...
                                             sum(insertInds), 1);
    end;
else
    disp('Computing maps from spatial basis weights');
    rePackedData=parcelFlag*data;
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