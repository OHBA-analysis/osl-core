function [simData, dipoleSignals,data] = simulate_MEG_signal(leadFields, ...
                                                        nSamples, ...
                                                        nChannels, ...
                                                        nTrials, ...
                                                        nDipoles, ...
                                                        fSample, ...
                                                        dipMeshInd, ...
                                                        dipoleOrientations, ...
                                                        signals, ...
                                                        SNR, ...
                                                        noiseLevel, ...
                                                        noiseBW, ...
                                                        roomNoiseCov, ...
                                                        roomNoiseSNR)
%SIMULATE_MEG_SIGNAL calculates a simulated signal in MEG sensors
%
% [SIMDATA, SIG] = SIMULATE_MEG_SINGAL(...) [all inputs required; not listed]
%    returns the signal in sensors, SIMDATA, based on a set of dipoles with
%    defined orientations and a set of leadfields for a brain mesh. The
%    signals input at each dipole location are returned in SIG. 
%
%    If more than one trial is specified, SIMDATA will hold each trial in a
%    separate cell. 
% 
%    White noise is added for each trial, determined by the supplied SNR or
%    a johnson shot noise amplitude and bandwidth
%    
%    Structured noise can be added if the covariance matrix is provided
%    (e.g. from measurements in an empty room) and an SNR
%
% See also: FT_DIPOLESIMULATION. 


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


%% Input processing
% check we have the right number of dipole locations
assert(isequal(length(dipMeshInd), nDipoles), ...
       [mfilename ':IncorrectDipMeshInd'], ...
       'Number of dipole indices on the mesh did not match nDipoles. \n');
   
% check lead-field size matches channels
assert(isequal(size(leadFields{1},1), nChannels), ...
       [mfilename ':InconsistentChannelsAndLFs'], ...
       'Number of channels did not match size of lead fields. \n');

% generate signals if not already existing
if isstruct(signals), 
    t     = (0:(nSamples - 1)) ./ fSample;
    phase = mod(signals.phase, 2*pi);
    s     = signals.amplitude .* cos(2*pi*signals.frequency.*t + phase);
    
    [dipoleSignals{1:nTrials}] = deal(repmat(s(:)', nDipoles, 1));
    
elseif nTrials==1 && ismatrix(signals), 
    dipoleSignals = {signals};
    
elseif iscell(signals),
    dipoleSignals = signals;
    
else %should not be here
    error([mfilename ':UnrecognisedSignalInput'], ...
          'Signal input not recognised. \n');
end%if

%% Simulate data
% calculate signals in sensors from each dipole
for iDipole = nDipoles:-1:1,
    dipoleMom = dipoleOrientations(iDipole, :)';
    
    %leadfield
    lf = leadFields{dipMeshInd(iDipole)};
    
    % multiply lead fields by dipole signals, weighted by orientation
    for iTrial = nTrials:-1:1,        
        data(iDipole).trial{iTrial} = zeros(nChannels, nSamples);
        for i = 1:3,
            data(iDipole).trial{iTrial} = data(iDipole).trial{iTrial} ...
                                   + lf(:,i) * ...
                                     (repmat(dipoleMom(i), 1, nSamples) .* ...
                                      dipoleSignals{iTrial}(iDipole, :));
        end%loop over orientations
    end%loop over trials
end%loop over dipoles

% Merge all N dipoles into one dataset
tic;
dipsum = data(1).trial;
if nDipoles > 1,
    % sum signals from each trial
    for iDipole = 2:nDipoles,
        for iTrial = 1:nTrials,
            newdata        = cell2mat(data(iDipole).trial(iTrial));
            dataprev       = cell2mat(dipsum(iTrial));
            dataraw        = newdata + dataprev;
            dipsum(iTrial) = mat2cell(dataraw, ...
                                      size(newdata,1), size(newdata,2));
        end%loop over trials
    end%loop over dipoles
end%if
toc;
simData       = data(1);
simData.trial = dipsum;

%% Add noise
for iTrial = 1:nTrials,
    dataraw   = cell2mat(simData.trial(iTrial));
    
    % calculate variance for noise specified as shot noise and/or SNR
    dataPower       = mean(dataraw(:).^2);
    if SNR, % prevent dividing by zero
        SNRnoiseVar = dataPower / SNR;
    else
        SNRnoiseVar = 0;
    end%if
    JohnsonNoiseVar = noiseLevel.^2 * noiseBW;
    noiseVar        = SNRnoiseVar + JohnsonNoiseVar;    
    
    % calculate white noise
    chanNoise = randn(size(dataraw)) .* sqrt(noiseVar);
    
    % calculate room noise
    % this only works for scanners where we have empty room data. The
    % amplitude of the room noise should be passed in as zero if it is not
    % applicable in this case
       
    % provide cholesky upper triangular factor of empty room noise if necessary
    if roomNoiseSNR && ~isempty(roomNoiseCov),
        [roomNoiseChol, isNotPosDef] = chol(roomNoiseCov);
        
        % check positive definite and size matches data
        if isNotPosDef,
            error([mfilename ':CovNotPosDef'], ...
                  ['Supplied noise covariance matrix must be ', ...
                   'positive definite. \n']);
        elseif ~isequal(size(roomNoiseCov, 1), size(dataraw, 1)),
            error([mfilename ':CovWrongSize'], ...
                  ['Supplied noise covariance matrix must match the ', ...
                   'number of simulated data channels. \n']);
        end%if
        
        roomNoise = MEGsim.randnorm(size(dataraw, 2), ...
                                    zeros(size(dataraw, 1), 1), ...
                                    roomNoiseChol);
                                
        % scale the covariance based on the SNR passed in
        roomNoisePower  = mean(roomNoise(:).^2);
        scaledRoomNoise = roomNoise .* ...
                          (dataPower / (roomNoiseSNR * roomNoisePower));
    else
        scaledRoomNoise = zeros(size(dataraw));
    end %if using room noise        
    
    % add noise
    dataraw = dataraw + chanNoise + scaledRoomNoise;
    
    % return to a cell
    simData.trial(iTrial) = mat2cell(dataraw, ...
                                     size(dataraw,1), ...
                                     size(dataraw,2));
end%loop over trials
end%simulate_MEG_signal
% [EOF]