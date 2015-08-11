function correlationMats = do_group_level_statistics(correlationMats, Settings)
%DO_GROUP_LEVEL_STATISTICS add group level inference 
% 
% CORRMATS = DO_GROUP_LEVEL_STATISTICS(CORRMATS, SETTINGS) adds group-level
%   inference fields to the CORRMATS structure. 
%
%   Adds fields groupEnvCorrelation_z, groupEnvPartialCorrelation_z and
%   groupEnvPartialCorrelationRegularized_z. 
%
%   Also adds structure falseDiscoveryRate, holding the standard z-stat
%   threshold over which network edge classifications will have a false
%   discovery rate of Settings.FDRalpha, and adjusted probabilities for
%   each network edge. Type `help ROInets.false_discovery_rate' for more
%   information. 
%
%   A fixed effects (z-test) or mixed effects (t-test) analysis on the mean 
%   z-stats over all subjects can be performed, set by the 
%   Settings.groupStatisticsMethod option. 
%
% See also: TTEST, ZTEST, ROInets.FALSE_DISCOVERY_RATE. 


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
%	$Revision: 231 $
%	$LastChangedDate: 2014-08-07 20:53:06 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 10-Apr-2014 13:33:21

if 1 == Settings.nSessions, 
    for iFreq = Settings.nFreqBands:-1:1,
        correlationMats{iFreq}.groupEnvCorrelation_z        = correlationMats{iFreq}.envCorrelation_z;
        correlationMats{iFreq}.groupEnvPartialCorrelation_z = correlationMats{iFreq}.envPartialCorrelation_z;
        if Settings.Regularize.do && isfield(correlationMats{iFreq}, 'envPartialCorrelationRegularized_z'),
            correlationMats{iFreq}.groupEnvPartialCorrelationRegularized_z = correlationMats{iFreq}.envPartialCorrelationRegularized_z;
        end%if
    end%for

else
    % take average correlations across subjects
    average = cellfun(@(C) ROInets.apply_function_to_structure(@take_average_correlations, C), ...
                      correlationMats, ...
                      'UniformOutput', false);

    switch lower(Settings.groupStatisticsMethod)
        case 'fixed-effects'
            % do z-test on mean of z-stats. This is a group-level fixed-effects
            % analysis, testing that the population mean of the z-stats at each node is
            % non-zero, under the population variance of 1. Note this is
            % well-described: the single-subject z-stats do follow a very good normal
            % distribution under the null, as tested with null data!
            % significance at 5% is z=1.96 and above
            % We don't need a t-test: we do t-tests when we approximate the population
            % s.d. with the sample s.d.
            % We know the population s.d. under the null here!
            sd_of_means = 1.0 ./ sqrt(Settings.nSessions);
            for iFreq = Settings.nFreqBands:-1:1,
                % this is equivalent to testing (mean(z) - mu) / s,
                % for mu=0 and s=1/sqrt(N) under H0: null distribution of sample means.

                correlationMats{iFreq}.groupEnvCorrelation_z        = ...
                    (average{iFreq}.envCorrelation_z        ./ sd_of_means);

                correlationMats{iFreq}.groupEnvPartialCorrelation_z = ...
                    (average{iFreq}.envPartialCorrelation_z ./ sd_of_means);

                if (Settings.Regularize.do && ...
                    isfield(average{iFreq}, 'envPartialCorrelationRegularized_z')),

                    correlationMats{iFreq}.groupEnvPartialCorrelationRegularized_z = ...
                        average{iFreq}.envPartialCorrelationRegularized_z            ...
                        ./ sd_of_means;
                end%if
            end%for

        case 'mixed-effects'
            % Mixed effects group-level stats: reduces power but allows for
            % mistakes/errors in modelling of 1st level stats - e.g. if scaling to
            % z-stats has buggered up, a mixed-effects analysis is still robust.
            % we use a t-test in this situation
            nullMean = 0;
            alpha    = 0.05;
            tail     = 'both';

            for iFreq = Settings.nFreqBands:-1:1,
                [~, p] = ttest(shiftdim(correlationMats{iFreq}.envCorrelation_z, 2), ...
                               nullMean, alpha, tail);
                signZ  = sign(average{iFreq}.envCorrelation_z);
                % convert p-value to two-tailed z-stat
                correlationMats{iFreq}.groupEnvCorrelation_z = ...
                             ROInets.p_to_z_two_tailed(shiftdim(p), signZ);

                [~, p] = ttest(shiftdim(correlationMats{iFreq}.envPartialCorrelation_z, 2), ...
                               nullMean, alpha, tail);
                signZ  = sign(average{iFreq}.envPartialCorrelation_z);
                correlationMats{iFreq}.groupEnvPartialCorrelation_z = ...
                             ROInets.p_to_z_two_tailed(shiftdim(p), signZ);

                if (Settings.Regularize.do && ...
                    isfield(correlationMats{iFreq}, 'envPartialCorrelationRegularized_z')),

                    [~, p] = ttest(shiftdim(correlationMats{iFreq}.envPartialCorrelationRegularized_z, 2), ...
                                   nullMean, alpha, tail);
                    signZ  = sign(average{iFreq}.envPartialCorrelationRegularized_z);
                    correlationMats{iFreq}.groupEnvPartialCorrelationRegularized_z = ...
                             ROInets.p_to_z_two_tailed(shiftdim(p), signZ);
                end%if
            end%for

        otherwise
            error([mfilename ':UnrecognisedGroupStats'], ...
                'Unrecognised group stats method: %s\n', ...
                Settings.groupStatisticsMethod);
    end%switch
end%if

% put group results through false discovery rate conversion to account for
% multiple network edge comparisons
for iFreq = Settings.nFreqBands:-1:1,
    % desired false discovery rate threshold
    correlationMats{iFreq}.falseDiscoveryRate.alpha = Settings.FDRalpha;

    % calculate equivalent z-threshold, and adjusted probability values on
    % each edge
    [~,                                                                           ...
     correlationMats{iFreq}.falseDiscoveryRate.zThreshold.groupEnvCorrelation_z,  ...
     correlationMats{iFreq}.falseDiscoveryRate.pAdjusted.groupEnvCorrelation_z] = ... 
 ROInets.false_discovery_rate(abs(correlationMats{iFreq}.groupEnvCorrelation_z),       ...
                              Settings.FDRalpha);

    [~,                                                                                  ...
     correlationMats{iFreq}.falseDiscoveryRate.zThreshold.groupEnvPartialCorrelation_z,  ...
     correlationMats{iFreq}.falseDiscoveryRate.pAdjusted.groupEnvPartialCorrelation_z] = ... 
 ROInets.false_discovery_rate(abs(correlationMats{iFreq}.groupEnvPartialCorrelation_z),       ...
                              Settings.FDRalpha);

    if (Settings.Regularize.do && ...
                    isfield(correlationMats{iFreq}, 'envPartialCorrelationRegularized_z')),
        [~,                                                                                             ...
         correlationMats{iFreq}.falseDiscoveryRate.zThreshold.groupEnvPartialCorrelationRegularized_z,  ...
         correlationMats{iFreq}.falseDiscoveryRate.pAdjusted.groupEnvPartialCorrelationRegularized_z] = ... 
     ROInets.false_discovery_rate(abs(correlationMats{iFreq}.groupEnvPartialCorrelationRegularized_z),       ...
                                  Settings.FDRalpha);
    end%if
end%for

end%do_group_level_statistics



%%% SUBFUNCTION %%%

function y = take_average_correlations(x)
if isnumeric(x), 
    y = nanmean(real(x), 3); 
else
    y = x; 
end%if
end%take_average_correlations
% [EOF]