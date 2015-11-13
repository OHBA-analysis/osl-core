function out = check_inputs(Settings)
%CHECK_INPUTS	Checks properties of Settings structure
%
% Settings = ROInets.check_inputs(Settings) checks inputs and establishes 
%   default settings


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
%	$Revision: 374 $
%	$LastChangedDate: 2015-01-12 17:30:23 +0000 (Mon, 12 Jan 2015) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 18:23:47

%% Define object
Inputs               = inputParser;
Inputs.CaseSensitive = false;
Inputs.FunctionName  = 'run_individual_network_analysis';
Inputs.StructExpand  = true;  % If true, can pass parameter-value pairs in a struct
Inputs.KeepUnmatched = false; % If true, accept unexpected inputs


%% validation functions
numericValidFcn     = @(x) (isnumeric(x) && isscalar(x)             && ...
                            ~isnan(x)    && ~isinf(x));

probabilityValidFcn = @(x) (numericValidFcn(x) && 0 <= x && x <= 1);

% check time range
validTimeRangeCheck = @(t) validateattributes(t, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nonnegative', 'nondecreasing'}, ...
                       Inputs.FunctionName, 'timeRange');
timeRangeFcn        = @(x) isempty(x)   ||                     ...
                           iscell(x)    ||                     ...
                           (isvector(x) && (2 == length(x)) && ...
                            (x(2) > x(1)) && all(x >= 0));
 
% check nifti files
isParcellationFcn   = @(x) (ischar(x)  && (exist(x, 'file')          ||   ...
                                           exist([x '.nii'], 'file') ||   ...
                                           exist([x '.nii.gz'], 'file'))) ...
                           || isnumeric(x) || islogical(x);
                                       
% check chosen orthog method
isProtocolFcn = @(c) ischar(validatestring(c,                                 ...
                                           {'none', 'symmetric',              ...
                                            'closest', 'pairwise'},           ...
                                           Inputs.FunctionName,               ...
                                           'leakageCorrectionMethod'));
                                        
% check chosen time-course method
isMethodFcn = @(c) ischar(validatestring(c,                            ...
                                        {'PCA', 'mean',                ...
                                         'spatialBasis', 'peakVoxel'}, ...
                                        Inputs.FunctionName,           ...
                                        'timecourseCreationMethod'));
                                    
% check group stats method
isGroupStatsMethodFcn = @(c) ischar(validatestring(c,                   ...
                                                   {'mixed-effects',    ...
                                                    'fixed-effects'},   ...
                                                   Inputs.FunctionName, ...
                                                   'groupStatisticsMethod'));
                                               
% check session names
isSessionNameFcn = @(c) ischar(c) || (iscell(c) && all(cellfun(@ischar, c)));

%% Inputs
Inputs.addParamValue('spatialBasisSet', [], isParcellationFcn);

Inputs.addParamValue('outputDirectory', [], @ischar);

Inputs.addParamValue('timeRange', [], timeRangeFcn);

RegularizeDefault = struct('do',        false, ...
                           'method',       [], ...
                           'path',         [], ...
                           'adaptivePath', [], ...
                           'Prior', struct('a', [], 'b', []));
Inputs.addParamValue('Regularize', RegularizeDefault, @isstruct);

protocolDefault = 'closest';
Inputs.addParamValue('leakageCorrectionMethod', protocolDefault, isProtocolFcn);

empiricalSamplesDefault = 0;
Inputs.addParamValue('nEmpiricalSamples', empiricalSamplesDefault, numericValidFcn);

ARmodelOrderDefault = 0;
Inputs.addParamValue('ARmodelOrder', ARmodelOrderDefault, numericValidFcn);

EnvParamsDefault = struct('windowLength', 2, ...
                          'overlap',      0.6,       ...
                          'useFilter',    false);
Inputs.addParamValue('EnvelopeParams', EnvParamsDefault, @isstruct);

freqBandDefault = {[]};
Inputs.addParamValue('frequencyBands', freqBandDefault, @iscell);

tcsMethodDefault = 'PCA';
Inputs.addParamValue('timecourseCreationMethod', tcsMethodDefault, isMethodFcn);

groupStatsMethodDefault = 'mixed-effects';
Inputs.addParamValue('groupStatisticsMethod', groupStatsMethodDefault, isGroupStatsMethodFcn);

FDRalphaDefault = 0.05;
Inputs.addParamValue('FDRalpha', FDRalphaDefault, probabilityValidFcn);

SaveCorrectedDef = struct('timeCourses',   false, ...
                          'envelopes',     true,  ...
                          'variances',     true,  ...
                          'ROIweightings', false);
Inputs.addParamValue('SaveCorrected', SaveCorrectedDef, @isstruct);

Inputs.addParamValue('sessionName', 'session1', isSessionNameFcn);
Inputs.addParamValue('gridStep',    8,          numericValidFcn);

%% Parse inputs
Inputs.parse(Settings);

% check any required inputs
if any(strcmpi('spatialBasisSet', Inputs.UsingDefaults)), 
    error([mfilename ':ParcelFileNotFound'], ...
          'You must provide a parcellation file or spatial basis set. \n');
end%if

% catch any settings not yet tailored for
if strcmpi('ica', Inputs.Results.spatialBasisSet),
    error('Performing the network analysis using the ICA results within the OIL structure is not yet possible. \n');
end%if

% check time range parameters
if ~isempty(Inputs.Results.timeRange),
    if iscell(Inputs.Results.timeRange),
        assert(isequal(length(Inputs.Results.timeRange), 1),                ... % only one session at a time for this function.
               [Inputs.FunctionName ':BadTimeRangeCellLength'],             ...
               ['%s: timeRange input must be empty, a time range as a ',    ...
                'two-component vector, or a cell array of two-component ',  ...
                'vectors. The cell array must have the same length as the', ...
                'number of subjects to do. \n'],                            ...
               Inputs.FunctionName);
        cellfun(validTimeRangeCheck, Inputs.Results.timeRange);
    else
        validTimeRangeCheck(Inputs.Results.timeRange);
    end%if
end%if

% check enveloping parameters
EnvelopeParams = Inputs.Results.EnvelopeParams;
numFieldNames = {'windowLength', 'overlap'};
for iff = 1:2,
    fieldName = numFieldNames{iff};
    if isfield(EnvelopeParams, fieldName),
        val = EnvelopeParams.(fieldName);
        assert(numericValidFcn(val) && val >= 0,            ...
               [mfilename ':IncorrectEnvField:' fieldName], ...
               ['EnvelopeParams.' fieldName ' must be a non-negative scalar. \n']);
    else
        EnvelopeParams.(fieldName) = EnvParamsDefault.(fieldName);
    end%if
end%for
if isfield(EnvelopeParams, 'useFilter'),
        assert(islogical(EnvelopeParams.useFilter),        ...
               [mfilename ':IncorrectEnvField:useFilter'], ...
               'EnvelopeParams.useFilter must be a logical value. \n');
else
    EnvelopeParams.useFilter = EnvParamsDefault.useFilter;
end%if
EnvelopeParams.useHanningWindow = false; % not an option at the moment due to a bug in osl

% check save request params
ff = fieldnames(SaveCorrectedDef);
SaveCorrected = Inputs.Results.SaveCorrected;
for iff = 1:length(ff),
    fieldName = ff{iff};
    if isfield(SaveCorrected, fieldName),
        assert(islogical(SaveCorrected.(fieldName)),                  ...
               [mfilename ':IncorrectSaveCorrectedField:' fieldName], ...
               ['SaveCorrected.' fieldName ' must be a logical value. \n']);
    else
        SaveCorrected.(fieldName) = SaveCorrectedDef.(fieldName);
    end%if
end%for


% check regularization parameters
Reg = Inputs.Results.Regularize;
assert(isfield(Reg, 'do') && islogical(Reg.do), ...
       [mfilename ':IncorrectRegularizationSpecification'], ...
       'The Regularize input must be a structure with a logical field ''do''. \n');
if Reg.do,
    if strcmpi(Inputs.Results.leakageCorrectionMethod, 'pairwise'),
        Inputs.Results.Regularize.do = false;
        warning([mfilename ':noRegularizationForPairwise'],                    ...
                ['It is not possible to compute regularized partial ',         ...
                 'correlation using the pairwise leakage correction method, ', ...
                 'yet. \n   The regularization has been turned off. \n']);
    else
        allowedRegMethods = {'Bayesian', 'Friedman'};
        assert((isfield(Reg, 'method')                            ...
                && any(strcmp(Reg.method, allowedRegMethods))),   ...
               [mfilename ':IncorrectRegularizationMethod'],      ...
               ['Regularize.method must be specified as either ', ...
                '''Bayesian'' or ''Friedman''. \n']);
           
           if strcmp(Reg.method, 'Bayesian'),
               isPriorFormatted = @(x) numericValidFcn(x) && x >= 0;
               isPriorSpecified = @(Reg) isstruct(Reg.Prior)             ...
                                  && all(isfield(Reg.Prior, {'a', 'b'})) ...
                                  && isPriorFormatted(Reg.Prior.a)       ...
                                  && isPriorFormatted(Reg.Prior.b);
               if isfield(Reg, 'Prior') && ~isempty(Reg.Prior),
                   assert(isPriorSpecified(Reg),                            ...
                         [mfilename, 'IncorrectRegularizationPrior'], ...
                         'Check the specification of Regularize.Prior. \n');
               else % use default Kerman neutral Gamma prior (2011)
                   fprintf(['  %s: Bayesian regularization selected. ',       ...
                            'Using default Kerman neutral Gamma hyperprior ', ...
                            'on regularization parameter. \n'], mfilename);
                   Reg.Prior = struct('a', 1/3, 'b', 0);
               end%if
                      
           else % Friedman
               assert(isfield(Reg, 'path') && ~isempty(Reg.path) && isnumeric(Reg.path), ...
                      [mfilename ':IncorrectRegularizationPath'],                        ...
                      'Please input one or more regularization parameters in Regularize.path. \n');
               assert(isfield(Reg, 'adaptivePath') && islogical(Reg.adaptivePath), ...
                      [mfilename ':IncorrectAdaptive'],                            ...
                      'Pleast specify a logical parameter in Regularize.adaptivePath. \n');  
               assert(all(Reg.path >= 0) && all(diff(Reg.path) > 0),      ...
                      [mfilename ':IncorrectRegularizationPathValues'],   ... 
                      ['Regression parameters must be non-negative and ', ...
                       'increasing in value. \n']);
           end%if
    end%if
end%if


% add any further useful options
Settings                     = Inputs.Results;
Settings.Regularize          = Reg;                                        % incorporate defaults
Settings.dateRun             = datestr(now, 31);
Settings.nSessions           = 1;
Settings.nFreqBands          = length(Inputs.Results.frequencyBands);
Settings.EnvelopeParams      = EnvelopeParams;                             % incorporate defaults
Settings.SaveCorrected       = SaveCorrected;                              % incorporate defaults



%% Output
out = Settings;
end%check_inputs
% [EOF]