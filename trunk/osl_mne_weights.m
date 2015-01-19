function [W, W_nonorm, lf] = osl_mne_weights(allLF, Cmu, Ce)
%OSL_MNE_WEIGHTS  computes MNE weights
%
%  Computes scalar weights for MNE source reconstruction, adjusting for
%  depth bias using sLORETA. 
%
% [W, W_NONORM, LF] = OSL_MNE_WEIGHTS(FULL_LF, DATACOV, NOISECOV) uses
%	lead fields mSensors x (nDipoles * 3) FULL_LF, the covariance of the
%	data over the sensors DATACOV and the covariance of room noise over the
%	sensors NOISECOV to estimate the weights vector used for calculation of
%	single-dipole magnitudes at each dipole location. 
%
%   Lead fields in FULL_LF must have three consecutive columns for each
%   dipole. Lead fields will be projected onto the single direction of
%   maximum variance for each dipole, independently. 
%
%   The function returns scalar weights with depth normalisation W in a
%   cell array of the same length as the number of dipoles, un-normalised
%   scalar weights W_NONORM and the projected lead fields LF. 
%
%   Data- and noise-covariance matrices should be pre-filtered for the
%   frequency band of interest. 	


%	References:
%	Dale, A.M. & Sereno, M.I.(1993) "Improved localization of cortical
%	activity by combining EEG and MEG with MRI cortical surface
%	reconstruction: a linear approach," J. Cogn. Neurosci 5, pp. 162--176. 
%
%	Pascual-Marqui, R.D. (2002) "Standardized low-resolution brain
%	electromagnetic tomography (sLORETA): technical details," Methods Find.
%	Exp. Clin. Pharmacol. 24 (Suppl D) pp. 5--12.
%
%   Hamalainen, M.S., Lin, F. & Mosher, J.C. (2010) "Anatomically and
%   functionally constrained minimum-norm estimates." In Hansen, P.C.,
%   Kringelbach, M.L. & Salmelin, R. (ed.) "MEG: an introduction to
%   methods," Oxford University Press, Oxford, UK, pp. 186--215. 

%	Copyright 2015 OHBA
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
%	Contact: giles.colclough@magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 19-Jan-2015 14:40:17

% dimensionality of inputs
mSensors = ROInets.rows(allLF);
D        = 3; % dimensions in space
nDipoles = ROInets.cols(allLF) ./ D;

% basic input size checking
assert(ROInets.rows(Ce) == mSensors && ROInets.cols(Ce) == mSensors, ...
       [mfilename ':InconsistentNoiseSensors'],                      ...
       'Rows of lead fields and size of noise covariance do not match. \n');
assert(ROInets.rows(Cmu) == mSensors && ROInets.cols(Cmu) == mSensors, ...
       [mfilename ':InconsistentDataSensors'],                         ...
       'Rows of lead fields and size of data covariance do not match. \n');
assert(~mod(nDipoles, 1),                       ...
       [mfilename ':InconsistentLFDimensions'], ...
       'Lead fields not provided in three dimensions. \n');

   
% regularization parameter
k = sum(eig(allLF * allLF.', Ce)) ./ (sum(eig(Cmu, Ce)) - mSensors); % sum(eig(B, A)) = trace(inv(A) * B);

% Full weights
Wfull = (allLF.') / (k * Ce + allLF * allLF.'); % A / B = A * inv(B)

% Find weights projected onto direction of maximal variance for each dipole
for iVox = nDipoles:-1:1,                                                  %#ok<BDSCI> - possible error caught above
    dipInd         = (D*iVox - (D-1)):(D*iVox);  % relevant indices for this dipole
    Ws             = Wfull(dipInd,:);            % single-dipole weights
    dipCov         = Ws * Cmu * Ws.';
    [ns,~]         = eigs(dipCov, [], 1);        % ns is direction of maximum source variance
    W_nonorm{iVox} = ns.' * Ws;                  % weights     for this dipole projected onto direction of maximum variance
    lf{iVox}       = allLF(:,dipInd) * ns;       % lead fields for this dipole projected onto direction of maximum variance
    lambda_s       = sqrt(W_nonorm{iVox} * (lf{iVox} * lf{iVox}.' ./ k + Ce) * W_nonorm{iVox}.'); % normalising constant for depth bias using sLORETA
    W{iVox}        = W_nonorm{iVox} ./ lambda_s; % normalised scalar weights vector for this dipole
end%loop over voxels

end%osl_mne_weights
% [EOF]