function [dipoleMag, neuralActivity] = LCMV_beamformer(Y, ...
                                                       L, ...
                                                       gridStep, ...
                                                       niiFileName, ...
                                                       saveNAI)
%LCMV_BEAMFORMER simple LCMV beamformer applied to demeaned, filtered data
%
% [DIPMAG, NAI] = LCMV_BEAMFORMER(Y, L)
%    beamforms demeaned, band-filtered data Y based on lead fields L for 
%    each voxel in the brain. L is a 1xnVoxels cell array. Each cell should
%    be an nSensorsx3 leadfield matrix. 
% 
%    The function returns DIPMAG, the dipole magnitudes estimated at each 
%    brain voxel, and NAI, the neural activity index over the voxels
%    inside the brain. 
%
% [...] = LCMV_BEAMFORMER(..., GRIDSTEP, FNAME, true) saves the neural 
%    activity index NAI in a .nii file specified by sting FNAME. 
%

% Reference: 
%   Van Veen, B.D.; Van Drongelen, W.; Yuchtman, M.; Suzuki, A., 
%   "Localization of brain electrical activity via linearly constrained 
%   minimum variance spatial filtering," Biomedical Engineering, IEEE 
%   Transactions on , vol.44, no.9, pp.867,880, Sept. 1997 
%   doi: 10.1109/10.623056


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


% process inputs
if nargin<5 || ~exist('saveNAI', 'var') || isempty(saveNAI),
    saveNAI = false;
end%if

fprintf('Beamforming data using LCMV. \n');

% find weights matrix by beamforming
C = cov(transpose(Y));
Cinv = pinv(C);

for iVoxel = length(L):-1:1,
    H         = L{iVoxel};
    M         = H' * Cinv * H;
    % find rank of M >=1
    rankM     = rank(M); if ~rankM, rankM = 1; end 
    W{iVoxel} = MEGsim.mwpinv(M, rankM) * H' * Cinv; % probably need to use mwpinv(M, 2) as M has rank 2.
    
    neuralActivity(iVoxel) = trace(MEGsim.mwpinv(M, rankM))  ...
                             ./ trace(MEGsim.mwpinv(H' * H, rankM));
end

% construct a lead field for each cartesian direction
for r = 3:-1:1,
    WTranspose{r} = transpose(cell2mat(cellfun(@(x) x(r,:)', ...
                              W, ...
                              'UniformOutput', false)));
end%loop over cartesian axes


fprintf('Finding magnitudes. \n');

% find dipole magnitudes
X         = cellfun(@(WT) WT*Y, WTranspose, 'UniformOutput', false);
catX      = cat(3,X{:});
dipoleMag = sqrt(sum(catX.^2,3));


% save power as .nii
if saveNAI,
    nii.quicksave(neuralActivity'./mean(neuralActivity), ...
                  niiFileName, ...
                  gridStep);
end%if
end%LCMV_beamformer