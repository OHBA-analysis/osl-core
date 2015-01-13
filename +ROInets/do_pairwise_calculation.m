function mats = do_pairwise_calculation(nodeData,        ...
                                        time,            ...
                                        EnvelopeParams)
% DO_PAIRWISE_CALCULATION pairwise source-leakage corrected network matrix
%
% MATS = DO_PAIRWISE_CALCULATION(DATA, TIME, ENVPARAMS) computes the node
%    correlation, envelope correlation and envelope partial correlations
%    for nodes corrected all-to-alll in a pairwise orthogonalisation
%    manner. 
%
%   DATA is nNodes x nTimeSamples, TIME is a vector of timepoints and
%   ENVPARAMS is a structure with fields:
%   - windowLength       : length of moving average window in s
%   - overlap            : fractional overlap of moving average window
%   - useHanningWindow   : use of a Hanning window function in the moving
%                          average block
%   - useFilter          : true/false: use moving average or more
%                          sophisticated downsampling / resampling filter. 
%    The passed frequency is 1.0/WINDOWLENGTH, resampled to a Nyquist 
%    frequency of twice this. 
%
%   MATS is a structure containing the Correlation, envCorrelation, and
%   envPartialCorrelation between nodes after the correction is applied. It
%   also contains the nSamples in the enveloped data, and a Regularization
%   field for tying into the osl_network_analysis pipeline. 
%
% extension of voxelwise orthogonalisation methods to ROI time-courses. 
% See Brookes et al 2012


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


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 373 $
%	$LastChangedDate: 2015-01-12 16:36:30 +0000 (Mon, 12 Jan 2015) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 18-Mar-2014 09:43:55

% do the calculation using pairwise orthogonalisation and partial
% correlation

nNodes = ROInets.rows(nodeData);

for iNode = nNodes:-1:1,
    fprintf('     Pairwise correlation analysis for ROI %d out of %d. \n', ...
            iNode, nNodes);
        
    otherNodes = ROInets.setdiff_pos_int(1:nNodes, iNode);
    % envelope now to save repeating inside the next loop
    iNodeEnv   = ROInets.envelope_data(nodeData(iNode,:), ...
                                       time,              ...
                                       EnvelopeParams);
    
    for ioindex = 1:length(otherNodes),
        oNode = otherNodes(ioindex);
        % remove the node timecourse from all others
        tmp = ROInets.householder_orthogonalise([nodeData(iNode,:)', ...
                                                 nodeData(oNode,:)']); % first vector is unaffected
        % orthNode will hold nodeData(oNode,:) pairwise orthogonalised to
        % iNode. It's transposed to a column vector. 
        orthNode = tmp(:,2);
        clear tmp
        
        % now remove these two from all others to get pcorr
        otherNodesSaveTwo = ROInets.setdiff_pos_int((1:nNodes), ...
                                                    [iNode,oNode]);
        
        dataForPCorr          = zeros(size(nodeData));
        dataForPCorr(iNode,:) = nodeData(iNode,:);
        dataForPCorr(oNode,:) = orthNode;
        
        % loop over all other nodes save these two
        for aoindex = 1:length(otherNodesSaveTwo),
            aoNode = otherNodesSaveTwo(aoindex);
            tmp    = ROInets.householder_orthogonalise([nodeData(iNode,:)', ...
                                                        orthNode, ...
                                                        nodeData(aoNode,:)'])'; % first two vectors are already orthogonal
            dataForPCorr(aoNode,:) = tmp(3,:);
        end%for
        
        % envelope
        envsForTwoNodes(otherNodes,:) = ROInets.envelope_data(dataForPCorr(otherNodes,:),    ...
                                                time,            ...
                                                EnvelopeParams);
        envsForTwoNodes(iNode,:) = iNodeEnv;
        
        % partial correlation
        covariance              = cov(envsForTwoNodes');
        partialCorrSingleRun    = ROInets.convert_precision_to_pcorr(...
                                     pinv(covariance));
        nodePCorr(iNode, oNode) = partialCorrSingleRun(iNode, oNode);
        
        % marginal correlation
        corrTmp = corrcov(covariance([iNode,oNode], [iNode,oNode])); % corr([envsForTwoNodes(iNode,:)', envsForTwoNodes(oNode,:)']);
        nodeCorr(iNode, oNode) = corrTmp(1,2);
        
        % raw correlation
        corrTmp = corr([nodeData(iNode,:)', orthNode]);
        rawCorr(iNode, oNode) = corrTmp(1,2);
    end%loop over other nodes
    
    nodeCorr(iNode,iNode)  = 1;
    nodePCorr(iNode,iNode) = 1;
    rawCorr(iNode,iNode)   = 1;
end%loop over nodes

mats.nSamples              = ROInets.cols(envsForTwoNodes); % this, unhelpfully, should be the number of samples in the enveloped data
mats.correlation           = rawCorr;
mats.envCorrelation        = nodeCorr;
mats.envPartialCorrelation = nodePCorr;
mats.Regularization.mean   = 0;
end%do_pairwise_calculation
% [EOF]