function [cleanD, jumpTimes] = osl_remove_jumps(D, varargin)
%OSL_REMOVE_JUMPS remove discontinuities from MEG raw signal
%
% [CLEAND, JUMPTIMES] = OSL_REMOVE_JUMPS(D) removes discontinuities from
%   SPM object D, by finding samples in the signal gradient which exceed
%   50 standard deviations of that gradient. Discontinuities are
%   interpolated by maintaining a smooth gradient. JUMPTIMES holds times of
%   discovered discontinuities to the nearest 0.1 s. CLEAND is the filtered
%   data, with additional events marking the jump locations. 
%
% [...] = OSL_REMOVE_JUMPS(D, 'PARAM', VALUE) takes additional information:
%   thresholdValue - set value for classification threshold [80]
%   thresholdType  - set type of threshold [std]:
%                        'abs': absolute change in gradient of signal 
%                        'std': change relative to std of gradient of signal
%                        'prctile': change relative to 95 percentile of signal
%   channels       - specify particular channels to filter [all MEG channels]
%   remove         - boolean flag to perform interpolation and
%                    identification of bad epochs in the data [true]
%   verbose        - boolean flag to report number of jumps found to
%                    standard out [true] 
%       
% This function removes "jumps" (discontinuities) from the EEG/MEG raw
% signal, based on a thresholding process, and filters the signal derivative
% over 20 timepoints.
% Such jumps occur with squid resetting and when acquisition is stopped
% with the "abort" button.
% This procedure is necessary before performing highpass filtering on the
% continuous data.
% Timestamps for the jumps are returned and at the same time recorded as
% event markers. Epochs containing such jumps should be rejected as they are
% affected by ringing from analogue filters in the recording system.
%   
%   See also SPM_EEG_REMOVE_JUMPS

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


%   $LastChangedBy$
%   $Revision$
%   $LastChangedDate$
%   Contact: giles.colclough@magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 27-Oct-2014 11:32:06

%% Input processing
% load in data, even if provided as a character string
D = spm_eeg_load(D);

% parse inputs
[threshVal, threshType, channels, Is] = assign_inputs(varargin{:});

if isempty(channels),
    channels = D.indchantype('MEGANY');
end%if

% create a new file
newFileName = fullfile(D.path, sprintf('j%s', fname(D)));
cleanD      = MEGsim.copy_meeg(D, newFileName);

% we will look for jumps in blocks of memory. 
[nBlocks, blockSize] = get_block_size(D, numel(channels));


%% filter blocks of channels
chanCount = 1;
for iBlock = 1:nBlocks,
    % load original data blockwise into workspace
    [Dtemp, chanCount, blockChanInds] = load_channel_block(D,         ...
                                                           chanCount, ...
                                                           channels,  ...
                                                           blockSize, ...
                                                           iBlock);
    
    jumps_fixed = false;
    
    % loop through channels within blocks
    for iChan = numel(blockChanInds):-1:1, % loop backwards to automatically initialise
        
        % find jumps in derivative
        data       = Dtemp(iChan, :, 1);
        dataDeriv  = diff(data);
        jumps      = find_jumps(dataDeriv, threshVal, threshType);
        
        % replace data
        if ~isempty(jumps) && Is.remove,
            
            % collapse jumps than are closer than 15 timepoints apart
            if numel(jumps) > 2, 
                jumps(find(diff(jumps) < 15) + 1) = []; 
            end%if

            data           = reconstruct_timecourse(data, dataDeriv, jumps);
            Dtemp(iChan,:) = data;
            jumps_fixed    = true;
        end%if
        
        % store jump onsets and filtered data
        if isequal(D.type, 'continuous')
            storedJumps{blockChanInds(iChan), 1}      = jumps;
        else
            storedJumps{blockChanInds(iChan), iBlock} = jumps;
        end%if
    end%for

    if jumps_fixed,
        % write filtered data blockwise in new data file
        if isequal(D.type, 'continuous')
            cleanD(channels(blockChanInds), :, 1) = Dtemp;
        else
            cleanD(channels, :, iBlock)           = Dtemp;
        end%if
    end%if
end%for

%% Insert artefact timepoints as event markers of type "artefact" 
% if remove, also label a bad epoch spanning 0.2 s.
[cleanD, jumpTimes] = label_jumps(cleanD, storedJumps, Is.remove);

%% Save new meeg object
cleanD = update_history(cleanD, D.path, D.fname, threshVal, ...
                        threshType, channels, Is.remove);
save(cleanD);
report_findings(fname(D), jumpTimes, Is.verbose);
end%osl_remove_jumps

%% Subfunctions
%--------------------------------------------------------------------------
function [thresh, type, channels, Is] = assign_inputs(varargin)
% input parser

% set up object
Inputs               = inputParser;
Inputs.CaseSensitive = false;
Inputs.FunctionName  = mfilename;
Inputs.StructExpand  = true;  % If true, can pass parameter-value pairs in a struct
Inputs.KeepUnmatched = false; % If true, accept unexpected inputs

% checking functions
numericValidFcn     = @(x) (isnumeric(x) &&  isscalar(x)  && ...
                            ~isnan(x)    &&  ~isinf(x));
channelValidFcn     = @(x)  isempty(x)   || (isnumeric(x) && ...
                            isvector(x)  && all(~mod(x,1)));
                        
% valid param-value options
Inputs.addParamValue('thresholdValue', 80, numericValidFcn);
Inputs.addParamValue('thresholdType', 'std', @(c) ischar(c));
Inputs.addParamValue('channels',       [], channelValidFcn);
Inputs.addParamValue('remove',         true, @(b) islogical(b));
Inputs.addParamValue('verbose',        true, @(b) islogical(b));

% run parsing
Inputs.parse(varargin{:});
validatestring(Inputs.Results.thresholdType, {'abs', 'std', 'prctile'}, ...
               mfilename, 'thresholdType');
           
thresh     = Inputs.Results.thresholdValue;
type       = Inputs.Results.thresholdType;
channels   = Inputs.Results.channels;
Is.remove  = Inputs.Results.remove;
Is.verbose = Inputs.Results.verbose;
end%assign_inputs

%--------------------------------------------------------------------------
function [nBlocks, blockSize] = get_block_size(D, nChannels)
%GET_BLOCK_SIZE find size of blocks of data to use

if isequal(D.type, 'continuous')
    % determine block size, depending on available memory
    if ispc
        % 2/3 of largest block of contiguous memory, for Windows platforms
        memsz = 2.0/3.0 * feature('memstats');
    else
        % 20 MB for all other platforms
        memsz = 20*1024*1024;
    end%if
    datasz    = nChannels * nsamples(D) * 8;  % datapoints x 8 bytes
    nBlocks   = ceil(datasz ./ memsz);
    blockSize = ceil(nChannels ./ nBlocks);
else
    nBlocks   = D.ntrials;
    blockSize = [];
end%if
end%get_block_size

%--------------------------------------------------------------------------
function [Dtemp, chanCount, blockChanInds] = load_channel_block(D, chanCount, channels, blockSize, iBlock)
%LOAD_CHANNEL_BLOCK loads a block of channels from D
if isequal(D.type, 'continuous'),
    blockChanInds = chanCount:(min(numel(channels), ...
                                   chanCount + blockSize - 1));
    Dtemp         = D(channels(blockChanInds), :, 1);
    chanCount     = chanCount + blockSize;

else
    blockChanInds = 1:length(channels);
    Dtemp         = D(channels, :, iBlock);
end
end%load_channel_block

%--------------------------------------------------------------------------
function jumps = find_jumps(dataDeriv, threshVal, threshType)
%FIND_JUMPS finds jumps in derivative of data under different thresholding
%  conditions
switch lower(threshType),
    case 'abs'
        threshold = threshVal;
        
    case 'std'
        threshold = threshVal * std(dataDeriv);
        
    case 'prctile'
        threshold = threshVal * prctile(abs(dataDeriv), 95);
        
    otherwise
        error([mfilename ':UnsupportedThresholdingMethod'], ...
              'Thresholding method %s not supported. \n',   ...
              threshType);
end%switch

jumps = find(abs(dataDeriv) > threshold);
end%find_jumps

%--------------------------------------------------------------------------
function data = reconstruct_timecourse(data, dataDeriv, jumps)
%RECONSTRUCT_TIMECOURSE fill in data around the jump

% need to wipe jump and some of the substantial ringing after
for iJump = 1:numel(jumps)
    % replace jump and timepoints -10 to +30 after
    replaceInd = jumps(iJump) + (-10:30);
    
    % calculate trend in signal from timepoints -30 to -5 after
    trendInd = jumps(iJump) + (-35:-5);
    
    % data might still be ringing when we link up - attempt to find mean
    % after jump and account for this
    newMeanInd = jumps(iJump) + (10:50);
    
    if replaceInd(1) < 1
        replaceInd = replaceInd + 1 - replaceInd(1);
        trendInd   = replaceInd;
    end%if
    
    if replaceInd(end) > length(dataDeriv)
        replaceInd = replaceInd - (replaceInd(end) - length(dataDeriv));
        trendInd   = replaceInd;
        newMeanInd = replaceInd(end) + 1; % this means no mean adjustment will happen
    end%if
    
    % interpolate around the jump
    dataDeriv(replaceInd) = mean(dataDeriv(trendInd));
    
    % account for mean afterwards
    remainingMeanShift         = mean(data(newMeanInd)) - ...
                                 data(replaceInd(end) + 1); % NB length(data) = length(dataDeriv) + 1, so this shouldn't break. 
    dataDeriv(replaceInd(end)) = dataDeriv(replaceInd(end)) - remainingMeanShift;
end%for

% reconstruct data
data = [0, cumsum(dataDeriv)] + data(1);

end%reconstruct_timecourse

%--------------------------------------------------------------------------
function [D, alljumps] = label_jumps(D, storedJumps, remove)
%LABEL_JUMPS find jump times in 0.1 second bins and label in events field
alljumps = cell(1, D.ntrials);

for iTrial = 1:D.ntrials
    % summarize jumps across channels into 0.1 s timebins, expressed in seconds
    alljumps{iTrial} = unique(ceil(cell2mat(storedJumps(:, iTrial).') ...
                                    / fsample(D) * 10))               ...
                        / 10.0;

    % find time of trial onset
    trialonset  = D.trialonset(iTrial);
    if iscell(trialonset)
        trialonset = trialonset{1};
    end%if
    
    if isempty(trialonset)
        trialonset = D.timeonset;
    end%if

    % extract events
    ev = events(D, iTrial);

    if iscell(ev)
        ev = ev{1};
    end%if

    nEvents = numel(ev);
    nJumps  = numel(alljumps{iTrial});
    % label jump events
    for iJump = 1:nJumps,
        ev(nEvents + iJump).type     = 'artefact';
        ev(nEvents + iJump).value    = 'jump';
        ev(nEvents + iJump).duration = [];
        ev(nEvents + iJump).time     = alljumps{iTrial}(iJump) + trialonset;
    end%for
    
    if remove,
        % mark bad epochs to cover ringing after the jump from hardware filters
        % run from 0.2 s before to 0.1 s after.
        for iJump = 1:nJumps,
            ev(nEvents + nJumps + iJump).type     = 'BadEpoch';
            ev(nEvents + nJumps + iJump).value    = 1;
            ev(nEvents + nJumps + iJump).duration = 0.3; %s
            ev(nEvents + nJumps + iJump).time     = alljumps{iTrial}(iJump) ...
                                                    - 0.2 + trialonset;
        end%for
    end%if
    if ~isempty(ev)
        [~, I] = sort([ev.time]);
        ev = ev(I);
        D = events(D, iTrial, ev);
    end%if
end%for

if numel(alljumps) == 1
    alljumps = alljumps{1};
end%if
end%label_jumps

%--------------------------------------------------------------------------
function D = update_history(D, oldPath, oldFname, threshVal, threshType, channels, remove)
%UPDATE_HISTORY updates history of altered D object
inputArgs = struct('D', fullfile(oldPath, oldFname), ...
                   'thresholdValue', threshVal,      ...
                   'thresholdType',  threshType,     ...
                   'channels',       channels,       ...
                   'remove',         remove);

D = history(D, mfilename, inputArgs);
end%update_history

%--------------------------------------------------------------------------
function [] = report_findings(dataFileName, jumpTimes, verbose)
%REPORT_FINDINGS
if ~verbose, return; end

nJumpsTotal = length(spm_vec(jumpTimes));

if 1 == nJumpsTotal, 
    gmEnding = ' ';
else
    gmEnding = 's';
end%if
fprintf('%s: Found %d jump%s in %s.\n', ...
        mfilename, nJumpsTotal, gmEnding, dataFileName); 
end%report_findings
% [EOF]