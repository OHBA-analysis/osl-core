function corrMats = convert_correlations_to_normal_variables(corrMats, ...
                                                             sigma,    ...
                                                             doRegularize)
%CONVERT_CORRELATIONS_TO_NORMAL_VARIABLES  converts correlations to z-scores
%
% CORRELATION_MATS = convert_correlations_to_normal_variables(CORRELATION_MATS, SIGMA, DO_REGULARIZE) 


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


%	$LastChangedBy: adambaker86@gmail.com $
%	$Revision: 261 $
%	$LastChangedDate: 2014-10-20 20:19:04 +0100 (Mon, 20 Oct 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Apr-2014 11:10:31

nullHypothesisCorrelation           = 0;

[~, corrMats.env_z]                 = ROInets.Fisher_r_to_z(corrMats.envCorrelation,   ...
                                                            corrMats.nSamples,         ...
                                                            nullHypothesisCorrelation, ...
                                                            sigma.z);
[~, corrMats.env_z_partial]         = ROInets.Fisher_r_to_z(corrMats.envPartialCorrelation, ...
                                                            corrMats.nSamples,              ...
                                                            nullHypothesisCorrelation,      ...
                                                            sigma.z_partial);

if doRegularize && isfield(corrMats, 'envPartialCorrelationRegularized'),
     % we use the unregularised width - best we can do. 
     % if we use the regularised width, then this is often compressed to
     % zero. 
    [~, corrMats.env_z_partial_reg] = ROInets.Fisher_r_to_z(corrMats.envPartialCorrelationRegularized, ...
                                                            corrMats.nSamples,                         ...
                                                            nullHypothesisCorrelation,                 ...
                                                            sigma.z_partial);
end%if
end%convert_correlations_to_normal_variables