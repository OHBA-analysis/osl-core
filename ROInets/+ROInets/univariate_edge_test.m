function [T, p, corrp] = univariate_edge_test(netmats, designMatrix, contrasts, standardise)
%UNIVARIATE_EDGE_TEST permutation test for group-level stats on networks
%
% [T, P, CORRP] = UNIVARIATE_EDGE_TEST(NETMATS, DESIGN, CONTRAST) performs
%    univariate testing for significance on each edge of a network matrix. 
%    Testing is performed using FSL's randomise, and uses 5000 permutations
%    of the group labels to perform nonparametric inference. 
%
%    Pass in NETMATS, which are symmetric network matrices with subjects in
%    the third dimension. The diagonals will be ignored. The DESIGN matrix
%    should have as many rows as subjects. The CONTRAST matrix should have
%    as many columns as the design matrix has columns. 
%
% [T, P, CORRP] = UNIVARIATE_EDGE_TEST(..., STANDARDISE) demeans and
%    variance normalises the design matrix, if TRUE. 
%
%   T provides single-edge T-stats; P the uncorrected p-Values; and CORRP the
%   FWE-corrected p-values. To use the weaker FDR correction, see
%   ROINETS.FALSE_DISCOVERY_RATE. 

%	Copyright 2015 OHBA, FMRIB
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


%	$LastChangedBy: GilesColclough $
%	$Revision: 763 $
%	$LastChangedDate: 2015-10-21 11:52:19 +0100 (Wed, 21 Oct 2015) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 25-Sep-2014 15:20:18


%% Input checking
[nNodes, checkMe, nSessions] = size(netmats);
assert(checkMe == nNodes, ...
       [mfilename ':NonsquareInput'], ...
       'Input netmats must be square, with subjects in the third dimension. \n');
   
[checkMe, nEVs] = size(designMatrix);
assert(checkMe == nSessions, ...
       [mfilename ':BadDesign'], ...
       'Design matrix must have as many rows as subjects. \n');
   
[nContrasts, checkMe] = size(contrasts);
assert(checkMe == nEVs, ...
       [mfilename ':BadContrasts'], ...
       'Contrasts must have as many columns as EVs in the design matrix. \n');

if nargin < 4 || ~exist('standardise', 'var'),
    standardise = false;
else
    assert(islogical(standardise), ...
           [mfilename ':BadStandardise'], ...
           'Standardise input must be a logical value. \n');
end%if
   
resultsDir = tempdir;

%% Construct design matrix
if standardise, 
    % demean and variance normalise
    X = bsxfun(@rdivide, bsxfun(@minus, designMatrix, mean(designMatrix)), ...
                         std(designMatrix));
else
    X = designMatrix;
end%if

% save out
designFile = fullfile(resultsDir, 'univariate_edge_test_design.mat');
ROInets.save_vest(X, designFile);
Cd = onCleanup(@() delete(designFile));

%% Construct contrasts
contrastFile = fullfile(resultsDir, 'univariate_edge_test_design.con');
ROInets.save_vest(contrasts, contrastFile);
Cc = onCleanup(@() delete(contrastFile));

%% Save out edges into nifti
inputNifti = fullfile(resultsDir, 'network_edges.nii.gz');
edges      = ROInets.get_edges(netmats); % note this assumes symmetry
for iS = ROInets.cols(edges):-1:1,
    formattedEdges(:,1,1,iS) = edges(:,iS);
end
save_avw(formattedEdges, inputNifti, 'f', [1 1 1 1]);
Ci = onCleanup(@() delete(inputNifti));

%% Run randomise
outputNifti = fullfile(resultsDir, 'univariate_edge_test');

% call to randomise
command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
                  inputNifti, outputNifti, designFile, contrastFile);
              
% submit to terminal
fprintf('%s\n', command);
[status, result] = system(command);
if status, 
    error([mfilename ':systemCallFailed'], ...
          'System command failed with message: \n   %s \n', result);
end%if 

Co = onCleanup(@() delete(outputNifti));

%% Retrieve results
for iCon = nContrasts:-1:1,
   TstatFile{iCon}     = [outputNifti '_tstat' num2str(iCon) '.nii.gz'];
   pFile{iCon}         = [outputNifti '_vox_p_tstat' num2str(iCon) '.nii.gz'];
   corrpFile{iCon}     = [outputNifti '_vox_corrp_tstat' num2str(iCon) '.nii.gz'];
   
   Ttmp(:,iCon)     = read_avw (TstatFile{iCon});
   ptmp(:,iCon)     = read_avw (pFile{iCon});
   corrptmp(:,iCon) = read_avw (corrpFile{iCon});
end%for

% tidy
for iCon = 1:nContrasts,
    delete(TstatFile{iCon});
    delete(pFile{iCon});
    delete(corrpFile{iCon});
end%for

% convert back to symmetric matrices
T = ROInets.unvectorize(Ttmp);
p = ROInets.unvectorize(1 - ptmp);
corrp = ROInets.unvectorize(1 - corrptmp);

end%univariate_edge_tests
% [EOF]