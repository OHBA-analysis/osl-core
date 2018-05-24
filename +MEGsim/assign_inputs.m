function out = assign_inputs(varargin)
%ASSIGN_INPUTS parses the inputs to parent function OSL_simulate_MEG_data


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
%	$Revision: 213 $
%	$LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45

%% Version check

% some people (Mark) still seem to have very old matlab ... newer versions
% are available on horus
if verLessThan('Matlab', '7.11'),
    warning([mfilename ':PleaseUpgrade'], ...
            ['Your version of matlab may not support functions like ', ...
             'ismatrix.m, which may produce errors. \n']);
end%if

%% Define object
Inputs               = inputParser;
Inputs.CaseSensitive = false;
Inputs.FunctionName  = 'osl_simulate_MEG_data';
Inputs.StructExpand  = false; % Can pass parameter-value pairs in a struct

%% validation functions
numericValidFcn     = @(x) (isnumeric(x) && isscalar(x)               && ...
                            ~isnan(x)    && ~isinf(x))                || ...
                           isempty(x);
                  
% check position is set of length-3 vectors
positionValidFcn    = @(x) ismatrix(x)  && isequal(size(x,2), 3); 

% check for MEEG object or filename
meegValidFcn        = @(x) (ischar(x)   && exist([x '.mat'], 'file')) || ...
                           isa(x, 'meeg');

% check signals are matrix or cell or structure
allowedSignalFields = {'amplitude'; 'frequency'; 'nTrials'; ...
                       'phase'; 'trialLength'};
                   
isSignalStruct      = @(x) isstruct(x)  && ...
                           all(strcmp(allowedSignalFields, ...
                                      fieldnames(orderfields(x))));
                                         
signalValidFcn      = @(x) (iscell(x)   && ...
                               (size(x,1) == 1 || size(x,2) ==1))     || ...
                           ismatrix(x)  || isSignalStruct(x);
                       
% check nifti files
isNiftiFcn          = @(x) isempty(x)   || (ischar(x)                 && ...
                                           (exist(x, 'file')          || ...
                                            exist([x '.nii'], 'file') || ...
                                            exist([x '.nii.gz'], 'file')));
                                        
% check covariance matrix
isCovarianceFcn     = @(x) isempty(x)   || (ismatrix(x)               && ...
                                            MEGsim.isposdef(x));
                            
                   
%% Inputs
% Mandatory inputs
addRequired(Inputs, 'dipolePositions',  positionValidFcn);
addRequired(Inputs, 'dipoleSignals',    signalValidFcn);
addRequired(Inputs, 'templateMEEGfile', meegValidFcn);

% Optional inputs

% Parameter-value pairs
fSampleDef = 250; % Hz
addParamValue(Inputs, 'fSample', fSampleDef, numericValidFcn);

dipoleOrientationsDef = [];
addParamValue(Inputs, 'dipoleOrientations', dipoleOrientationsDef, ...
              @(x) positionValidFcn(x) || isempty(x));
          
resDef = 2; % mm
addParamValue(Inputs, 'spatialResolution',  resDef,   numericValidFcn);

SNRDef = 0; 
noiseLevelDef = 0; % sensor units / sqrt(Hz)
noiseBWDef    = 0; % Hz
roomNoiseCovDef = []; 
roomNoiseSNR    = 0;

addParamValue(Inputs, 'whiteSignalToNoiseRatio', SNRDef,   numericValidFcn);
addParamValue(Inputs, 'noiseLevel',         noiseLevelDef, numericValidFcn);
addParamValue(Inputs, 'noiseBandWidth',     noiseBWDef,    numericValidFcn);
addParamValue(Inputs, 'emptyRoomNoiseCovariance', ...
              roomNoiseCovDef, isCovarianceFcn);
addParamValue(Inputs, 'structuredSignalToNoiseRatio', ...
              roomNoiseSNR,    numericValidFcn);


bfDef = false;
addParamValue(Inputs, 'runBeamformer',        bfDef,      @(x) islogical(x));

debugDef = false;
addParamValue(Inputs, 'debugMode',            debugDef,   @(x) islogical(x));

allChanDef = false;
addParamValue(Inputs, 'matchTemplateExactly', allChanDef, @(x) islogical(x));

rndSeedDef = sum(100*clock); % standard practice for initialisation
addParamValue(Inputs, 'randomNumberSeed',   rndSeedDef, numericValidFcn);

fileNameDef = [tempname(pwd) '_meg'];

while 1,
    if exist(fileNameDef, 'file'), 
        fileNameDef = [fileNameDef(1:end-1), ...
                       num2str(str2double(fileNameDef(end))+1)];
    else%valid filename
        break;
    end
end%while
addParamValue(Inputs, 'fileName', fileNameDef, @(x) ischar(x));

blankDef = []; % inherit from base MEEG file
addParamValue(Inputs, 'fiducials',         blankDef, ...
              @(x) isstruct(x) || isempty(x));
addParamValue(Inputs, 'channels',          blankDef, ...
              @(x) isstruct(x) || isempty(x));
addParamValue(Inputs, 'sensors',           blankDef, ...
              @(x) isstruct(x) || isempty(x));
addParamValue(Inputs, 'trials',            blankDef, ...
              @(x) isstruct(x) || isempty(x));
addParamValue(Inputs, 'beamformerResults', blankDef, ...
              @(x) isstruct(x) || isempty(x)); % BF structure from running source recon
addParamValue(Inputs, 'sourceReconSaveFile', blankDef, ...
              @(x) ischar(x) || isempty(x));   % file to save BF results to
addParamValue(Inputs, 'modalities', blankDef, ...
              @(x) ischar(x) || isempty(x));   % 
          
addParamValue(Inputs, 'structuralFile', '', isNiftiFcn)
          
%% Parse inputs
parse(Inputs, varargin{:});

% deduce number of trials from signal input
if isstruct(Inputs.Results.dipoleSignals),
    nTrials = Inputs.Results.dipoleSignals.nTrials;
    
elseif iscell(Inputs.Results.dipoleSignals),
    nTrials = length(Inputs.Results.dipoleSignals);
    
else % single matrix input
    nTrials = 1;
    
end%if

% deduce number of dipoles from position input
nDipoles = size(Inputs.Results.dipolePositions, 1);

% check no confusion over dipole entry
if isequal(size(Inputs.Results.dipolePositions, 2), nDipoles),
    warning([mfilename ':uncertainDipolePositions'], ...
            ['3 Dipoles specified. ', ...
             'Ensure each dipole position forms a row. \n']);
end%if

% check signals match nDipoles and no confusion over input 
if iscell(Inputs.Results.dipoleSignals),
    signalnRows    = cellfun(@(C) size(C, 1), Inputs.Results.dipoleSignals);
    signalnSamples = cellfun(@(C) size(C, 2), Inputs.Results.dipoleSignals);
    nSamples       = signalnSamples(1); % check for consistency later
    
    isCorrectSize  = all((signalnRows(:)    == 1)         | ...
                         (signalnRows(:)    == nDipoles));
    isSizeDubious  = all((signalnSamples(:) == nDipoles)  & ...
                         (signalnRows(:)    == nDipoles));
    isSamplesConst = all( signalnSamples    == nSamples);
    
elseif isstruct(Inputs.Results.dipoleSignals),
    % structure input -> signals handled later
    isCorrectSize  = true;
    isSizeDubious  = false;
    isSamplesConst = true; 
    nSamples       = round(Inputs.Results.fSample * ...
                           Inputs.Results.dipoleSignals.trialLength); 
        
elseif ismatrix(Inputs.Results.dipoleSignals),
    [signalnRows, nSamples] = size(Inputs.Results.dipoleSignals);
    isCorrectSize  = (signalnRows(:) == 1) | (signalnRows == nDipoles);
    isSizeDubious  = isequal(nSamples, signalnRows);
    isSamplesConst = true;
    
else
    error([mfilename ':unrecognisedSignalInput'], ...
           'Input for dipoleSignals should be struct, matrix or cell. \n');
end%if

if ~isCorrectSize,
    error([mfilename ':IncorrectNumberDipSignals'], ...
          ['The number of columns in each cell of dipoleSignals must be ', ...
           '1 or equal to the number of columns in dipolePositions. \n']);
end%if
if isSizeDubious,
    warning([mfilename ':uncertainDipoleSignals'], ...
            ['Length of signals in each cell match the number of dipoles. ', ...
            'Ensure each dipole signal forms a row. \n']);
end%if
if ~isSamplesConst,
    error([mfilename ':Inconsistent sample lengths'], ...
          ['The number of samples (rows) in each cell of dipoleSignals ', ...
           'must be the same. \n']);
end%if

% check orientations match nDipoles
if (~isempty(Inputs.Results.dipoleOrientations) ...
    && ~isequal(size(Inputs.Results.dipoleOrientations, 1), nDipoles)),
    error([mfilename ':IncorrectNumberDipOrientations'], ...
          ['The number of columns in dipoleOrientations and ', ...
           'dipolePositions must be the same. \n']);
end%if

% check trials input matches nTrials
if ~isempty(Inputs.Results.trials), 
    assert(isequal(length(Inputs.Results.trials), nTrials), ...
           [mfilename ':IncorrectNumberTrialEntries'], ...
           ['The number of fields in the input trial structure must ', ...
            'match the number of trials as specified by the signals entry. \n']);
end%if

% check noise bandwidth supplied if noise level is given
if ( Inputs.Results.noiseLevel && ~Inputs.Results.noiseBandWidth) || ...
   (~Inputs.Results.noiseLevel &&  Inputs.Results.noiseBandWidth),
    error([mfilename ':NoiseLevelAndBWMismatch'], ...
          ['When specifying an absolute noise level, please supply both ', ...
           'the noise level in sensor units per root Hz and the noise ', ...
           'bandwidth in Hz. \n']);
end%if

% check room noise covariance both supplied if using
if ( Inputs.Results.structuredSignalToNoiseRatio && ...
     isempty(Inputs.Results.emptyRoomNoiseCovariance)) || ...
   (~Inputs.Results.structuredSignalToNoiseRatio && ...
    ~isempty(Inputs.Results.emptyRoomNoiseCovariance)),
    error([mfilename ':RoomNoiseMismatch'], ...
          ['When specifying structured noise, please supply both ', ...
           'the noise covariance matrix and a related ', ...
           'signal-to-noise ratio  \n']);
end%if 

% check that room noise is supplied in a size which matches the data
if Inputs.Results.structuredSignalToNoiseRatio && ...
   ~isempty(Inputs.Results.emptyRoomNoiseCovariance),
    if isempty(Inputs.Results.channels),
        Dtemplate = spm_eeg_load(Inputs.Results.templateMEEGfile);
        nChannels = length(MEGsim.megchannels(Dtemplate));
    else
        nChannels = length(Inputs.Results.channels);
    end%if
    assert(isequal(nChannels, ...
                   size(Inputs.Results.emptyRoomNoiseCovariance, 1)), ...
           [mfilename ':RoomNoiseCovarianceAndTemplateMismatch'], ...
           ['The supplied noise covariance does not match the ', ...
            'number of channels to be simulated. \n']);
end%if using structured noise    

%% Allocate inputs to variables
out.nDipoles            = nDipoles;
out.nTrials             = nTrials;
out.nSamples            = nSamples;
out.fSample             = Inputs.Results.fSample;
out.dipolePositions     = Inputs.Results.dipolePositions;
out.dipoleOrientations  = Inputs.Results.dipoleOrientations;
out.signals             = Inputs.Results.dipoleSignals;
out.templateMEEGfile    = Inputs.Results.templateMEEGfile;
out.spatialRes          = Inputs.Results.spatialResolution;
out.SNR                 = Inputs.Results.whiteSignalToNoiseRatio;
out.noiseLevel          = Inputs.Results.noiseLevel;
out.noiseBW             = Inputs.Results.noiseBandWidth;
out.emptyRoomCovariance = Inputs.Results.emptyRoomNoiseCovariance;
out.roomNoiseSNR        = Inputs.Results.structuredSignalToNoiseRatio;
out.channels            = Inputs.Results.channels;
out.sensors             = Inputs.Results.sensors;
out.trials              = Inputs.Results.trials;
out.fiducials           = Inputs.Results.fiducials;
out.doBeamformer        = Inputs.Results.runBeamformer;
out.saveFileName        = Inputs.Results.fileName;
out.structuralFile      = Inputs.Results.structuralFile; 
out.beamformerResults   = Inputs.Results.beamformerResults;
out.reconResFileName    = Inputs.Results.sourceReconSaveFile;
out.DEBUG               = Inputs.Results.debugMode;
out.matchAllChannels    = Inputs.Results.matchTemplateExactly;
out.rndSeed             = Inputs.Results.randomNumberSeed;
out.modalities             = Inputs.Results.modalities;
end%assign_inputs
% [EOF]