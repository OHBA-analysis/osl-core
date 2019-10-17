function D = osl_erase_trigger_artefact(D, triggerChan, triggerTime, deleteWindow, trialInds, doPlot)
%OSL_ERASE_TRIGGER_ARTEFACT zeros out data around OHBA's 2015 trigger artefact
%
% Use this function to clean your images before publication, and to make it
% clear in any subsequent analyses if you've forgotten to exclude this
% window. Use it after any filtering processes you intend to do. 
% Alternatively, run your whole analyses, and zero out the artefact from
% your final TF/ERF plots.
%
% D = OSL_ERASE_TRIGGER_ARTEFACT(D, TRIGGERCHAN) zeros a window of 40ms 
%  around t = 0s  every trial. TRIGGERCHAN is the string giving the label
%  of the trigger channel
%
%  The function will clean the data, over-writing D. It will set the
%  zeroed-out data as artefact events, flagging the bad samples.
%
% D = OSL_ERASE_TRIGGER_ARTEFACT(D, TRIGGERCHAN, TRIGGERTIME, DELETEWINDOW, TRIALINDS)
%  optionally pass in a set of times to use, TRIGGERTIME, in seconds. 
%  Use one value (the same for each trial), or one value per trial. 
%  You can set the window size DELETEWINDOW in ms for zeroing (default 
%  40ms). You can also choose a subset of trials, specified by TRIALINDS.
%
% D = OSL_ERASE_TRIGGER_ARTEFACT(D, TRIGGERCHAN, [], ...) will attempt to
%  auto-detect the times to zero using data in the trigger channel given by
%  label TRIGGERCHAN. (e.g. 'STI101'). Check your data after using this
%  option: success is not guaranteed!
%
% D = OSL_ERASE_TRIGGER_ARTEFACT(..., DOPLOT) will provide rude
%  epoch-averaged plots to check the result.
%
% You should ensure that you exclude a short window centred on the trigger
% from all further analyses, so that these zeros are not included. 
%
% How to descibe this analysis in your methods / supplementary methods:
%   We have removed data associated with a trigger-locked artefact from our
%   figure for clarity. We ignored a window of X ms centred on the
%   trigger in all analyses. 



%   Copyright 2015 OHBA
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


%   $LastChangedBy: GilesColclough $
%   $Revision: 763 $
%   $LastChangedDate: 2015-10-21 11:52:19 +0100 (Wed, 21 Oct 2015) $
%   Contact: giles.colclough@magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 01-Dec-2015 14:30:16

%% Input checking
% D should be an epoched, sensor-space, Neuromag object
D = spm_eeg_load(D); 
assert(~strcmpi(D.type, 'continuous'), ...
       [mfilename ':NotEpochedD'], ...
       'Input spm object should be epoched data. \n');
assert(any(strcmp('MEGMAG', D.chantype)) && any(strcmp('MEGPLANAR', D.chantype)), ...
       [mfilename ':NotNeuromag'], ...
       'Input spm object should be a sensor-space Neuromag recording. \n');
   
% plot result?
if nargin < 6 || ~exist('doPlot', 'var'),
    doPlot = false;
else
    doPlot = logical(doPlot);
end%if

% select trials to use
if nargin < 5 || ~exist('trialInds', 'var') || isempty(trialInds),
    trialInds = 1:D.ntrials;
else
    assert(all(ismember(trialInds, 1:D.ntrials)), ...
           [mfilename ':BadTrialSelection'],       ...
           'Selected trials did not match number of trials in object. \n');
end%if

nTrialsToUse = length(trialInds);

% set size of window to zero out
if nargin < 4 || ~exist('deleteWindow', 'var') || isempty(deleteWindow),
    deleteWindow = 40; %ms
else
    assert(isnumeric(deleteWindow) && all([1,1] == size(deleteWindow)), ...
           [mfilename ':BadWindowSpec'],                                ...
           'deleteWindow should be a time in ms. \n');
end%if
% check that the window matches length of trial
trialLength = max(D.time) - min(D.time); % in s
assert(deleteWindow./1000 < trialLength, ...
       [mfilename ':WindowTooLong'],     ...
       'Window to zero is longer than trial length. \n');

% trigger chan should be the name of one of the channels in D 
assert(ischar(triggerChan) && any(strcmp(triggerChan, D.chanlabels)), ...
       [mfilename ':BadChanInput'],                                   ...
       'TRIGGERCHAN should be the name of a channel in D. \n');
   
if nargin < 3 || ~exist('triggerTime', 'var'),
    triggerTime = 0;
end%if
if isempty(triggerTime),
    autoDetect  = true;
    triggerTime = 0;
else
    autoDetect = false;
end%if
if all(1 == size(triggerTime,1) && 1 == size(triggerTime,2)),
    triggerTime = repmat(triggerTime, nTrialsToUse, 1);
end%if
assert(all(triggerTime < max(D.time) & triggerTime > min(D.time)), ...
       [mfilename ':BadTriggerTimes'],                             ...
       'Trigger times should be within time range of trial. \n');

    




%% Preliminaries
% pull out helpful channels
triggerInd = find(strcmp(triggerChan, chanlabels(D)));
megInd     = find(strncmpi('MEG', chantype(D), 3)); % spm's function routinely misses channels

% have a look at data before cleanup
if doPlot,
    hf = figure('Color', 'w', 'Name', 'Trigger artefact before cleanup');
    plot_problem(hf, D, triggerInd, trialInds);
end%if




%% Zero relevant data
if autoDetect, 
    zeroInds = auto_detect_zero_inds(D, triggerInd, deleteWindow, trialInds);
else
    zeroInds = get_zero_inds(D, triggerTime, deleteWindow);
end%if

% we might have differing zero lengths in each trial
for iTrial = 1:nTrialsToUse,
    trial  = trialInds(iTrial);
    toZero = zeroInds{iTrial};
    D(megInd, toZero, trial) = zeros(length(megInd), length(toZero));
end%for

D = set_bad(D, zeroInds, trialInds);
                                   
                                   
                                   
%% Plot cleaned data and save result before returning
D.save;

if doPlot,
    hf = figure('Color', 'w', 'Name', 'Trigger artefact after cleanup');
    plot_problem(hf, D, triggerInd, trialInds);
end%if



end%osl_remove_trigger_artefact




function zeroInds = get_zero_inds(D, triggerTime, deleteWindow)
%GET_ZERO_INDS parse information to produce indices to zero

% number of zeros either side of trigger time
padding = round(deleteWindow./1000 * D.fsample ./ 2); 

for iTrial = length(triggerTime):-1:1,
    trigInd          = nearest(D.time, triggerTime(iTrial));
    zeroInds{iTrial} = trigInd-padding:trigInd+padding;
end%for

end%get_zero_inds

function zeroInds = auto_detect_zero_inds(D, triggerInd, deleteWindow, trialInds)
%AUTO_DETECT_ZERO_INDS finds indices to zero pased on trigger channel

% number of zeros either side of trigger time
padding = round(deleteWindow./1000 * D.fsample ./ 2); 

for iTrial = length(trialInds):-1:1,
    trigger = D(triggerInd,:,trialInds(iTrial));
    
    % assume the baseline in a trigger channel is its mode
    triggerOn = trigger ~= mode(trigger);
    
    diffOn    = [0 diff(triggerOn)];
    onSamples = find(1 == diffOn);
    for iS = 1:length(onSamples),
        triggerOn(onSamples-padding:onSamples+padding) = 1;
    end%for
    offSamples = find(-1 == diffOn);
    for iS = 1:length(offSamples),
        triggerOn(offSamples-padding:offSamples+padding) = 1;
    end%for
    zeroInds{iTrial} = find(triggerOn);
end%for

end%get_zero_inds

function D = set_bad(D, zeroInds, trialInds)
%SET_BAD marks artefact chunks

for iTrial = 1:length(trialInds),
    trial = trialInds(iTrial);
    % get old events
    Events = D.events(trial);
    
    BadEvents = struct([]);
    if ~isempty(zeroInds{iTrial}),
        BadEvents(1).type     = 'artefact_OHBA_trigger';
        BadEvents(1).value    = 'all';
        BadEvents(1).time     = D.time(zeroInds{iTrial}(1));
        BadEvents(1).duration = (zeroInds{iTrial}(end) - zeroInds{iTrial}(1) + 1) ...
                                     / D.fsample;
        BadEvents(1).offset   = 0;
    end%if
    
    % Concatenate new and old events
    if ~isempty(BadEvents)
        Events = [Events{:}; BadEvents(:)];
    end
    % Save new events with previous
    D = events(D,trial,Events);
end%for
end%set_bad

% make a quick function to create epoch averages
function average = epoch_average(D, chanInds, trialInds)
if strcmpi(D.type, 'continuous'),
    average = 0;
    for iTrial = 1:length(trialInds),
        trialSamples = trialInds{iTrial};
        average      = average + D(chanInds, trialSamples, :);
    end%for
    average = average ./ length(trialInds);
else
    average = mean(D(chanInds,:,trialInds),3);
end%if
end%epoch_average

% make a quick function to create plots
function [] = plot_problem(hf, D, triggerInd, trialInds)

% we've already checked for neuromag data so this is a safe assumption
magInds    = find(strcmp('MEGMAG',    D.chantype));
planarInds = find(strcmp('MEGPLANAR', D.chantype));

figure(hf);
subplot(3,1,1);
plot(D.time, epoch_average(D, magInds, trialInds)'); % D.time, 
% xlim([-0.4 0.6])
ylabel('fT');
subplot(3,1,2);
plot(D.time, epoch_average(D, planarInds, trialInds)'); % D.time, 
% xlim([-0.4 0.6])
ylabel('fT/mm');
subplot(3,1,3);
plot(D.time, epoch_average(D, triggerInd, trialInds)'); % D.time, 
% xlim([-0.4 0.6])
ylabel('uV');
xlabel('time / s');
drawnow;
end%plot_problem

% [EOF]