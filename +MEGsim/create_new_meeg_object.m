function Dout = create_new_meeg_object(Dtemplate, filename, P)
%CREATE_NEW_MEEG_OBJECT creates new meeg object by inheritance
%
% D = CREATE_NEW_MEEG_OBJECT(TEMPLATE, FILENAME, P) creates new MEEG object
%     D with FILENAME, using based on TEMPLATE and using parameters in P. 
%     
%     If P.trialData does not exist, or is empty, D holds blank data.
%     Otherwise, D holds data determined by fieldtrip data structure
%     P.TRIALDATA. 
%
% P must contain:
%   nTrials    - number of trials in the created object
%   nSamples   - number of time samples in the created object
%   nChannels  - number of channels in the created object
%   fSample    - sampling frequency in the created object
%
% P may contain:
%   trialData - 1xnTrials cell array. Each cell contains an nChannels x
%               nSamples MEG data matrix. If this field is not specified, a
%               blank MEEG object is created. 
%
%   simulatedChannelIndices - a list of indices specifying the channels in
%                             the template which the supplied data will
%                             fill. Default behaviour is to strip out all
%                             other channels
%
%   matchAllChannels        - if TRUE, all channels in the template are
%                             retained, and the supplied data replace the
%                             channels specified by simulatedChannelIndices
%
%   noCheckSens             - if TRUE, performs only a simple check on the
%                             D object. This is useful if an MEG dataset
%                             has used the EEG box for artefact channels,
%                             but has no information on the EEG sensors. 
%
% P may also be used to replace the event information, fiducials or sensors
% structures. 
%   sensors   - replaces the D.sensors('MEG') field
%   fiducials - replaced the D.fiducials field
%   trials    - fields 'label', 'trialonset', 'events' can be used to
%               replace the event information. 
%
% See also: MEEG/CLONE, SPM_EEG_MONTAGE, OSL_CHANGE_SPM_EEG_DATA,
% OSL_CONCAT_SPM_EEG_CHANS. 

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
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45

%% input check
if ~ischar(filename), 
    error([mfilename ':NonCharNewFileName'], ...
          'new file name should be a character string. \n');
end%if

timeOnsetDef = 0;

% check to see if only a subset of the template's channels are used
if isfield(P, 'simulatedChannelIndices') && ~isempty(P.simulatedChannelIndices)
    simChanInd  = P.simulatedChannelIndices(:);
    
    if ~isfield(P, 'matchAllChannels') || ~P.matchAllChannels,
        % nChannels will be reduced to number of simChanInd
        P.nChannels = length(simChanInd);
    else
        % channels will stay same as template, but we fill only a subset
        % with data.
    end%if
    
    assert(~islogical(simChanInd), ...
           [mfilename ':LogicalChanInds'], ...
           'Can''t use logical vectors to index channels in MEEGs. \n');
else
    simChanInd  = ':'; % use all channels
end%if

%% copy across data, or start afresh
% ensure that there is nothing already residing at filename
MEGsim.check_filename_is_clear(filename)
    
% If new data specified
isNewDataSpecified = isfield(P, 'trialData') && ~isempty(P.trialData);
if isNewDataSpecified,
    
    fprintf('Creating MEEG object with new data. \n');

    if (~isfield(P, 'matchAllChannels') || ~P.matchAllChannels), 
        D = create_new_D_from_scratch(Dtemplate, filename, P, ...
                                      simChanInd, timeOnsetDef);

    else % match all channels
        D = create_new_D_from_clone(Dtemplate, filename, P, ...
                                      simChanInd, timeOnsetDef);
    end%if match all channels

% if not specified, create blank datafile. 
else
  fprintf('Copying across template with blank data. \n');
  % clone. Consider resetting history
  D = clone(Dtemplate, filename, [P.nChannels, P.nSamples, P.nTrials]);
  
  D = fsample(D, P.fSample);
  D = timeonset(D, timeOnsetDef);
end%if

%% set channels
if isfield(P, 'channels') && ~isempty(P.channels),
    error(['Can''t change channel info atm - not tested. \n ', ...
           'Use a template of the scanner format you want to replicate. \n']);
    
    assert((P.nChannels == length(P.channels)), ...
           [mfilename ':ChannelMismatch'], ...
           ['Mismatch between number of channels and ', ...
            'length of channel structure. \n']);                           %#ok<UNRCH>
        
    D = badchannels(D, arrayfun(@(S) S.bad, P.channels, ...
                                'UniformOutput', false));
    D = units(D, arrayfun(@(S) S.units, P.channels, ...
                                'UniformOutput', false));
    D = chanlabels(D, (1:P.nChannels)', arrayfun(@(S) S.label, P.channels, ...
                                                 'UniformOutput', false));
    D = chantype(D, (1:P.nChannels)', arrayfun(@(S) S.type, P.channels, ...
                                               'UniformOutput', false));
                                           
    xCoor = cat(2, P.channels.X_plot2D);
    yCoor = cat(2, P.channels.Y_plot2D);
    D     = coor2D(D, (1:P.nChannels)', cat(1, xCoor, yCoor));        
end

%% set trials and events info
if isfield(P, 'trials') && ~isempty(P.trials),
    warning([mfilename ':EventSettingNotTested'], ...
            'Not tested with setting event info. Use at your own risk. \n');
    
    if isfield(P.trials, 'label'),
        cond = {P.trials.label};
        D    = conditions(D, [], cat(1, cond{:}));
    end%if    
    if isfield(P.trials, 'trialonset'),
        to = {P.trials.trialonset};
        D  = trialonset(D, ':', cat(1,to{:}));
    end%if
    if isfield(P.trials, 'events'),
        D = events(D, ':', {P.trials.events});
    end%if
else
    % Should already be taken care of during clone/copy. 
    % Checkmeeg will fill in the rest. 
end%if    

%% copy inv 
if isfield(Dtemplate, 'inv'),
    D.inv = Dtemplate.inv;
    D.save;
end

%% Set sensors
% set properties of MEG sensors. This should correspond to the output of
% ft_datatype_sens. The 'sensors' structure is also commonly called 'grad'. 
if isfield(P, 'sensors') && ~isempty(P.sensors),
    warning([mfilename ':SensorsSet'], ...
            'Not tested with new sensor properties. \n');
    D = sensors(D, 'MEG', P.sensors);
    
    % set duplicated fields
    if isfield(D, 'inv'),
        D.inv{D.val}.datareg.sensors = P.sensors;
    end%if
end%if 
if isempty(sensors(D, 'MEG')) && isNewDataSpecified,
    
    % consruct a new sensor structure using just the channels from the sim
    simulatedLabels = chanlabels(D, ':');
    grad            = sensors(Dtemplate, 'MEG');
    
    for iLabel = length(grad.label):-1:1, 
        matchingSensInd(iLabel) = any(strcmpi(grad.label{iLabel}, ...
                                              simulatedLabels));
    end
    matchingSensInd = logical(matchingSensInd);
    
    newgrad          = grad;
    newgrad.chantype = grad.chantype(matchingSensInd);
    newgrad.chanpos  = grad.chanpos(matchingSensInd,:);
    newgrad.chanori  = grad.chanori(matchingSensInd,:);
    newgrad.chanunit = grad.chanunit(matchingSensInd);
    newgrad.label    = grad.label(matchingSensInd);
    newgrad.tra      = grad.tra(matchingSensInd, :);
    
    % below taken straight from spm_eeg_ft2spm.m as way of setting sensors
    % on new object
    D = sensors(D, 'MEG', ft_convert_units(newgrad, 'mm'));
    
    S = [];
    S.task = 'project3D';
    S.modality = 'MEG';
    S.updatehistory = 0;
    S.D = D;

    D = spm_eeg_prep(S);
    
%     % note to use all sensors in template, even irrelevant ones, we would
%     % do:
%     D = sensors(D, 'MEG', sensors(Dtemplate, 'MEG'));                    error(); % error added in case uncommented by mistake
end%if

%% Set fiducials
if isfield(P, 'fiducials') && ~isempty(P.fiducials),
    warning([mfilename ':FiducialsSet'], ...
            'Not tested with new fiducial properties. \n');
    D = fiducials(D, P.fiducials);
    % set duplicated fields
    if isfield(D, 'inv'),
        D.inv{D.val}.datareg.fid_eeg = P.fiducials;
    end%if
end%if
if isempty(fiducials(D)) && isNewDataSpecified,
    D = fiducials(D, fiducials(Dtemplate));
end%if

%% Check object is ok
if isfield(P, 'noCheckSens') && P.noCheckSens, % e.g. helpful if have used EEG box to record artefact channels. EEG channels with no EEG sens information will cause errors in checkmeeg. 
    [Dout,res] = check(D);
else
    [Dout,res] = check(D, 'sensfid');
end

if ~res, 
    error([mfilename ':FailedMEEGCheck'], ...
          ['MEEG object is not internally consistent. ', ...
           'Hopefully spm has given you some more helpful errors... \n']);
end%if

% check data type
if ~strcmpi(type(Dout), 'continuous'), 
    warning([mfilename, ':notContinuousData'], ...
            ['This function has mostly been tested for continous data ', ...
             'generation. \n  ', ...
             'Please report any bugs when generating trial data. \n']);
end%if

res = save(Dout);
% if ~res, 
%     error([mfilename ':FailedMEEGSave'], ...
%           ['MEEG object has failed to save. ', ...
%            'Hopefully spm has given you some more helpful errors... \n']);
% end%if
end%create_new_meeg_object






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function D = create_new_D_from_scratch(Dtemplate, filename, P, ...
    simChanInd, timeOnsetDef)
%CREATE_NEW_D_FROM_SCRATCH creates MEEG object and fills with data

if ~iscell(P.trialData) && isnumeric(P.trialData) % 1 trial
    P.trialData = {P.trialData};
elseif ~iscell(P.trialData),
    error([mfilename ':IncorrectDataInput'], ...
          'Expecting a cell or numeric array of data in P.trialData. \n');
end%if


% check trial consistency
assert((length(P.trialData) == P.nTrials), ...
       [mfilename ':FTTrialMismatch'], ...
       ['Mismatch between number of trials in input and ', ...
        'number of trials in data. \n']);

% input new data
ftData.trial = P.trialData;

% use the same time vector for each trial
ftData.fsample = P.fSample;
time           = (timeOnsetDef:(P.nSamples - 1)) ./ P.fSample;
for iTrial = 1:P.nTrials,
    ftData.time{iTrial} = time;
end%for

% we use only those channels that were simulated
ftData.label = chanlabels(Dtemplate, simChanInd);

% check channel consistency

assert((P.nChannels == length(ftData.label)), ...
       [mfilename ':FTChannelMismatch'], ...
       ['Mismatch between number of channels in input and ', ...
        'number of channels in template. \n']);
    
assert((P.nChannels == size(ftData.trial{1}, 1)), ...
       [mfilename ':FTChannelMismatch2'], ...
       ['Mismatch between number of channels in input and ', ...
        'number of channels in data. \n']);

% possibly needed for neuromag support - see L140 of spm_eeg_ft2spm.m
templateGrad = sensors(Dtemplate, 'MEG');
hdr.Fs =  P.fSample;
hdr.grad = templateGrad;
hdr.label = templateGrad.label;

if ft_senstype(hdr, 'neuromag'),
    % ft_chantype fails to correctly identify channel types. Let's help
    % it out by providing it with the fields it expects
    hdr.orig = get_orig_struct_for_neuromag_labelling(hdr, Dtemplate);
    
    % include the hdr for the spm_eeg_ft2spm function
    ftData.hdr = hdr;
end%if neuromag

% create new meeg object. Should handle internal consistency
D = spm_eeg_ft2spm(ftData, filename);

% Seems to be an error with neuromag. Catch mis-labeling of channels
if strcmpi(modality(D), 'Other'),
    warning([mfilename ':ModalityOtherEvent'], ...
            ['Some spm functions appear to be failing to assign ', ...
             'correct channel labels. This may cause problems in oil ', ...
             'with spm_eeg_montage. \n The cause may be a neuromag ', ...
             'template - this is not well tested atm. ', ...
             'Try using ctf templates. \n']);
    D = chantype(D, (1:P.nChannels)', chantype(Dtemplate, simChanInd));
    D = units(D, (1:P.nChannels)', units(Dtemplate, simChanInd));
    D = coor2D(D, (1:P.nChannels)', coor2D(Dtemplate, simChanInd));
end%if
end%create_new_D_from_scratch
%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function D = create_new_D_from_clone(Dtemplate, filename, P, ...
                                     simChanInd, timeOnsetDef)
%CREATE_NEW_D_FROM_CLONE clones template and replaces data
if ~iscell(P.trialData) && isnumeric(P.trialData) % 1 trial
    P.trialData = {P.trialData};
end%if

% check trial consistency
assert((length(P.trialData) == P.nTrials), ...
       [mfilename ':FTTrialMismatch'], ...
       ['Mismatch between number of trials in input and ', ...
        'number of trials in data. \n']);

% check channel consistency
assert((length(simChanInd) == size(P.trialData{1}, 1)), ...
       [mfilename ':FTChannelMismatch'], ...
       ['Mismatch between number of channels simulated and ', ...
        'number of channels in data. \n']);

% clone sensor, channel and fiducial info
D = clone(Dtemplate, filename, [P.nChannels, P.nSamples, P.nTrials], 2); % clone sensor and fiducial info

D(simChanInd, :, :) = cat(3, P.trialData{:});
D = fsample(D, P.fSample);
D = timeonset(D, timeOnsetDef);

% to avoid errors in lead fields, we are aiming to use every channel in
% the template
assert(isequal(nchannels(Dtemplate), nchannels(D)), ...
       [mfilename ':ChannelMismatchBetweenTemplateAndNewObj'], ...
       ['Unless you use all the channels which where in the template, ', ...
        'I cannot guarantee that lead field calculations are accurate. \n']);

% Seems to be an error with neuromag. Catch mis-labeling of channels
if strcmpi(modality(D), 'Other'),
    warning([mfilename ':ModalityOtherEvent'], ...
            ['Some spm functions appear to be failing to assign ', ...
             'correct channel labels. This may cause problems in oil ', ...
             'with spm_eeg_montage. \n The cause may be a neuromag ', ...
             'template - this is not well tested atm. ', ...
             'Try using ctf templates. \n']);
    D = chantype(D, (1:P.nChannels)', chantype(Dtemplate, simChanInd));
    D = units(D, (1:P.nChannels)', units(Dtemplate, simChanInd));
    D = coor2D(D, (1:P.nChannels)', coor2D(Dtemplate, simChanInd));
end%if
% check we have sensor info for all channels
Dgrad = sensors(D, 'MEG');
assert(all(ismember(lower(chanlabels(D, MEGsim.megchannels(D))), ...
                    lower(Dgrad.label))), ...
       [mfilename ':NotEnoughSensorInfo'], ...
       'There is not sufficient sensor info to cover all channels. \n');
end%create_new_D_from_clone





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function orig = get_orig_struct_for_neuromag_labelling(hdr, Dtemplate)
%GET_ORIG_STRUCT_FOR_NEUROMAG_LABELLING
% ft_chantype fails to correctly identify channel types. Let's help
% it out by providing it with the fields it expects
%
% See also: FT_CHANTYPE (L91). 

orig.channames.KI = zeros(length(hdr.label),1); % holds channel kind, 1=meg
orig.chaninfo.TY  = zeros(length(hdr.label),1); % holds coil type (0=magnetometer, 1=planar gradiometer)

gradChanTypes = hdr.grad.chantype;

% grad chantypes can often be just 'unknown' in neuromag -
% this is unhelpful. We will provide a fix
if ~any(strncmpi('MEG', gradChanTypes, 3))
    templateChanTypes = chantype(Dtemplate, ':');
    for iLabel = length(gradChanTypes):-1:1,
        if (strcmpi('unknown', gradChanTypes(iLabel)) || ...
            strcmpi('other', gradChanTypes(iLabel))),
                
            matchingInd = strcmpi(hdr.label{iLabel}, ...
                                  chanlabels(Dtemplate));
                                       
            gradChanTypes(iLabel) = templateChanTypes(matchingInd);
        end%if
    end%loop over labels
end%if

% fill in labelling fields
orig.channames.KI(strncmpi('MEG', gradChanTypes, 3))   = 1;
orig.chaninfo.TY( strcmpi('MEGPLANAR', gradChanTypes)) = 1;
end%get_orig_struct_for_neuromag_labelling
% [EOF]