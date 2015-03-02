function [] = osl_example_ROI_network_analysis()
% Example ROI network analysis

% load in an oil structure after source recon
% Look at your data before you proceed with this function!
oilMat = fullfile('/path/to/oil/matrix', ...
                  'oil.mat');
fprintf('Loading oil matrix %s\n', oilMat);
tmp = load(oilMat, 'oil');
oil = tmp.oil;
clear tmp;

% choose a binary ROI map. Take care that the resolution of the nifti file
% matches that of the source reconstruction.
parcellationDirectory = '/path/to/parcellation/';
parcelFile = fullfile(parcellationDirectory, ...
                      'parcellationFile.nii.gz');
                  
% choose a results directory
outDir = '/path/to/results/';

% setup the ROI network settings
% type `help osl_network_analysis' for more detailed information
oil.ROInetworks = struct();
oil.ROInetworks.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
oil.ROInetworks.subjectsToDo             = [2,8,14,17,20,23,26,29] ;       % vector of indices indicating which of the source reconstruction sessions to analyse
oil.ROInetworks.nParallelWorkers         = 0;                              % analyse several sessions in parallel, if non-zero. Specifies number of matlab workers to use. Caution: these tasks are very memory-heavy! You might slow down because you max out on RAM. 
oil.ROInetworks.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
oil.ROInetworks.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
oil.ROInetworks.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
oil.ROInetworks.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
oil.ROInetworks.leakageCorrectionMethod  = 'symmetric';                    % choose from 'symmetric', 'pairwise' or 'none'. 
oil.ROInetworks.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
oil.ROInetworks.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
oil.ROInetworks.envelopeWindowLength     = 2; % s                          % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
oil.ROInetworks.envelopeWithFilter       = true;                           % use a more sophisticated filter than a sliding window average
oil.ROInetworks.frequencyBands           = {[8 13], [13 30], []};          % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
oil.ROInetworks.timecourseCreationMethod = 'spatialBasis';                 % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
oil.ROInetworks.outputDirectory          = outDir;                         % Set a directory for the results output
oil.ROInetworks.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
oil.ROInetworks.FDRalpha                 = 0.05;                           % false determination rate significance threshold


% run the ROI network analysis
[oil, networkMatrices] = osl_network_analysis(oil);

% save the network matrices in nifti files for easier viewing
% I tend to use a resampled parcellation file for this. It's prettier. 
prettyParcelFile =  fullfile(parcellationDirectory, ...
                             'upsampledParcellationFile.nii.gz');
prettyGridStep   = 2; % mm
imageSaveDir     = fullfile(oil.ROInetworks.outputDirectory, 'images', filesep);
ROInets.make_directory(imageSaveDir);
imageNameRoot    = fullfile(imageSaveDir, 'myNetworkMatrices');
BAND_NAMES       = {'Alpha', 'Beta', 'Broadband'};
save_network_mats_as_nifti(networkMatrices,  ...
                           prettyParcelFile, ...
                           prettyGridStep,   ...
                           imageNameRoot,    ...
                           BAND_NAMES,       ...
                           oil.ROInetworks.doRegularize);
end%osl_example_network_analysis









function [ ] = save_network_mats_as_nifti(netMats, ...
                                          parcelFileName, ...
                                          gridStep,       ...
                                          imageNameRoot,  ...
                                          BAND_NAMES,     ...
                                          doRegularize)
%SAVE_NETWORK_MATS_AS_NIFTI	Saves network matrices in nifti files
%
% [] = SAVE_NETWORK_MATS_AS_NIFTI(NETMATS, PARCELFILE, GRIDSTEP, ...
%                                 IMAGENAMEROOT, BAND_NAMES, DOREGULARIZE)
%   saves the z-statistic maps as nifti files for easy viewing. 
%   The z-statistic matrices in NETMATS are saved into nifti files, with
%   file name roots (including full path) specified by IMAGENAMEROOT. 
%   The parcellation file used for display is PARCELFILE, a nifti file
%   containing binary ROI allocations, with as many ROIs as the size of the
%   network matrices in NETMATS. GRIDSTEP is the resolution of this parcel
%   file. 
%   BAND_NAMES is a cell array of strings holding names for each frequency
%   band which forms a cell element of NETMATS. 
%   DOREGULARIZE indicates if the regularized partial correlation matrices
%   are present in NETMATS. 
%	
%   See also: OSL_NETWORK_ANALYSIS. 


%	Copyright 2014 OHBA
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
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 18-Mar-2014 11:25:01


if ~isempty(parcelFileName) && ischar(parcelFileName), 
    parcelMasks = nii_quickread(parcelFileName, gridStep);
else
    error([mfilename ':MaskReadError'], ...
          'Unable to read ROI mask file %s. \n', ...
          parcelFileName);
end%if

% use a binary parcel mask to plot the network matrices
parcelMasks = logical(parcelMasks);

nFreqs = length(netMats);
assert(isequal(nFreqs, length(BAND_NAMES)), ...
       [mfilename ':FreqBandMismatch'], ...
       ['Length of frequency band names does not match length of ', ...
        'network matrix cell. \n']);
    
workbenchDir = '/path/to/workbench/';

for iFreq = 1:nFreqs,
    orthEnvCorrSaveFileName  = [imageNameRoot '_' BAND_NAMES{iFreq} ...
                                '_envelope_correlations_z'];
    orthEnvPCorrSaveFileName = [imageNameRoot '_' BAND_NAMES{iFreq} ...
                                '_envelope_partial_correlations_z'];
    envRegSaveFileName       = [imageNameRoot '_' BAND_NAMES{iFreq} ...
                                '_envelope_partial_correlations_z_lasso'];
    
    ROInets.nii_parcel_quicksave(netMats{iFreq}.groupEnvCorrelation_z,        ...
                                 parcelMasks, [], orthEnvCorrSaveFileName,        ...
                                 gridStep);
    ROInets.nii_parcel_quicksave(netMats{iFreq}.groupEnvPartialCorrelation_z, ...
                                 parcelMasks, [], orthEnvPCorrSaveFileName,       ...
                                 gridStep);
                             
    osl_render4D([orthEnvCorrSaveFileName '.nii.gz'],  ...
                 orthEnvCorrSaveFileName, workbenchDir);
    osl_render4D([orthEnvPCorrSaveFileName '.nii.gz'], ...
                 orthEnvPCorrSaveFileName, workbenchDir);
    
    
    if (doRegularize && ...
        isfield(netMats{iFreq}, 'envPartialCorrelationRegularized_z')),
        ROInets.nii_parcel_quicksave(netMats{iFreq}.groupEnvPartialCorrelationRegularized_z, ...
                                     parcelMasks, [], envRegSaveFileName,                        ...
                                     gridStep);
        osl_render4D([envRegSaveFileName '.nii.gz'], ...
                     envRegSaveFileName, workbenchDir);
    end%if
end%for
end%osl_save_network_mats_as_nifti