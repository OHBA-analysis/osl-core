function sigma = find_empirical_H0_distribution_width(nH0Iter, nNodes, nTimeSamples, Settings, RegularizationResults, ARmodel, Fs, Filter, EnvelopeParams)
%FIND_EMPIRICAL_H0_DISTRIBUTON_WIDTH
%
% generation of empirical H0 dataset by using fitted AR model 
% and calculating distribution of null correlations.

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
%	$Revision: 263 $
%	$LastChangedDate: 2014-10-23 11:30:39 +0100 (Thu, 23 Oct 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 23:39:39



if nH0Iter, 
    fprintf(' Generating empirical null correlation distribution. \n');
    % set regularization to be applied to empirical data
    if Settings.Regularize.do,
        if strcmpi(Settings.Regularize.method, 'Friedman'),
            empirical_rho = RegularizationResults.mean;
        else % Bayesian
            % subsample from posterior of regularization parameters. If
            % nH0Iter is greater than the number of posterior samples,
            % choose with replacement. 
            nPosteriorSamples = length(RegularizationResults.posteriorRho);
            if nH0Iter <= nPosteriorSamples,
                replacement = false;
            else 
                replacement = true;
            end%if
            empirical_rho = RegularizationResults.posteriorRho(...
                                randsample(nPosteriorSamples,      ...
                                           nH0Iter, replacement));
        end%if
    else
        empirical_rho = 0;
    end%if

    sigma.nH0Iter  = nH0Iter;
    sigma.nNodes   = nNodes;
    sigma.nSamples = nTimeSamples;
    [sigma.z,         ...
     sigma.z_partial, ...
     sigma.z_partial_reg] = ROInets.find_empirical_dev(nH0Iter,        ...
                                                       sigma.nNodes,   ...
                                                       sigma.nSamples, ...
                                                       ARmodel,        ...
                                                       Settings.Regularize.do, ...
                                                       empirical_rho,  ... 
                                                       Fs,             ...
                                                       Filter,         ...
                                                       EnvelopeParams);
else
    % we will assume normality
    [sigma.z,         ...
     sigma.z_partial, ...
     sigma.z_partial_reg] = deal([]);
end%if
end%find_empirical_H0_distribution_width