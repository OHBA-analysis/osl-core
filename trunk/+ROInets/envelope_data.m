function [envelopedData, t_ds, newFs] = envelope_data(dipoleMag, t, varargin)
%ENVELOPE_DATA applies Hilbert envelope to data, without normalisation
%
% ENV = ENVELOPE_DATA(DIPOLEMAG, TIME) takes dipoles magnitudes in a set
%    of voxels over TIME, applies a hilbert envelope and downsamples using
%    a moving-average window length of 2 seconds, with no overlap. 
%    Dipole magnitudes are expected to occupy rows of DIPOLEMAG: (nVoxels x
%    nSamples). 
%
% ENV = ENVELOPE_DATA(..., 'windowLength', WINDOWLENGTH) uses 
%    WINDOWLENGTH seconds when downsampling using a moving average. 
%
% ENV = ENVELOPE_DATA(..., 'overlap', OVERLAP) uses a
%    fractional overlap OVERLAP between windows.
%
% ENV = ENVELOPE_DATA(..., 'useHanningWindow', true)
%    applies a Hanning envelope to the moving average window
%
% ENV = ENVELOPE_DATA(..., 'windowLength', WINDOWLENGTH, 'useFilter', true) 
%    uses a more sophisticated downsamping operation, with low-pass 
%    filtering and care over edge-effects. 
%    The passed frequency is 1.0/WINDOWLENGTH, resampled to a Nyquist 
%    frequency of twice this. 
%
% ENV = ENVELOPE_DAT(..., 'takeLogs', true) returns the logarithm of the
%    power envelopes, computed as 2*log(env).
%
% ENV = ENVELOPE_DATA(..., 'verbose', true) increases the volume of text
%    written to stdout. 
%
% [ENV, DS_TIME] = ENVELOPE_DATA(...) returns the downsample time vector
%    DS_TIME of same length as ENV.


%	Copyright 2013-4 OHBA
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
%	$Revision: 367 $
%	$LastChangedDate: 2014-12-13 19:04:35 +0000 (Sat, 13 Dec 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45

% inputs
[windowLength, overlap, useHanningWindow, useFilter, takeLogs, verbose] = ...
    parse_inputs(varargin{:});

% Demean
signal = ROInets.demean(dipoleMag, 2);
clear dipoleMag

% Hilbert Envelope with nested function
HE = hilbert_envelope(verbose);
clear signal
    
if verbose,
    fprintf('%s: downsampling. \n', mfilename);
end%if

% downsample and low-pass filter
if ~useFilter,
    % do osl-style averaging
    [envelopedData, t_ds, newFs] = moving_average_downsample(t,            ...
                                                             overlap,      ...
                                                             windowLength, ...
                                                             useHanningWindow);
    
else
    % use a better filter, with control on edge effects. Be aware that this
    % is different to published approaches for finding band-limited power.
    newMaxFreq            = 1.0 / windowLength;
    % use nested function to minimize memory movement
    [envelopedData, t_ds, newFs] = filter_and_downsample(t, newMaxFreq);
    
    % still getting some issues on the last sample. Let's throw it out. 
    envelopedData(:,end) = [];
    t_ds(end)            = [];
    
end%if useFilter

% convert to logarithm of power
if takeLogs,
    envelopedData = 2 .* log(real(envelopedData) + eps(min(abs(envelopedData(:))))); % prevent log(0).
end%if

if verbose, 
    fprintf('%s: complete\n', mfilename); 
end%if
%%% END of Main Function %%%







%%%%  SUBFUNCTIONS  %%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HE = hilbert_envelope(verbose)
%HILBERT_ENVELOPE computes Hilbert envelope of signal
%
% H = HILBERT_ENVELOPE(S) computes the Hilbert envelope, H, of
%   signal S. Each columns of S is a sample in time, rows contain
%   different signals. H is returned in the same format.
if verbose,
    fprintf('%s: finding Hilbert transform. \n', mfilename);
end%if

nSamples = ROInets.cols(signal);
nfft     = 2 ^ nextpow2(nSamples + 1); % zero-padding improves accuracy and speed
for iVox = ROInets.rows(signal):-1:1,
    HE(iVox,:) = abs(transpose(hilbert(transpose(signal(iVox,:)), nfft)));
    % free up memory
    signal(iVox,:) = [];
end%for
HE       = HE(:, 1:nSamples);
end%hilbert_envelope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [envelopedDataDS, t_ds, newFs] = filter_and_downsample(t, newMaxFreq)
%FILTER_AND_DOWNSAMPLE downsamples signal Matlab's decimate function
%
% [SIGDS, TDS] = FILTER_AND_DOWNSAMPLE(SIG, T, NEWFS)
%   downsamples signal SIG to new sampling rate NEWFS. 
%   Each column of SIG is a
%   sample in time, rows of SIG hold different signals. T is a vector
%   indexing the time at each sample of SIG. A constant sampling rate is
%   assumed. 
%
%   The function returns the downsampled signal, SIGDS, and a downsampled
%   time vector, TDS.
%
% See also: DECIMATE. 


% use matlab's decimate function, which lowpass filters and eliminates edge
% effects

Fs            = 1.0 / median(diff(t));
nVoxels       = ROInets.rows(HE);
nyquistFactor = 2.0 / 0.8;
newFs         = nyquistFactor * newMaxFreq;
DsFactor      = floor(Fs / newFs);
newFs         = Fs / DsFactor;
factorList    = get_factors(DsFactor);
filterType    = 'iir';%'fir'; % use fir filter for long time-series, and for higher downsample factors. Be aware that no phase shift correction is applied

% downsample in factors lower than 13.
for iFactor = 1:length(factorList),
    for iVoxel = nVoxels:-1:1, % loop backwards to initialise correctly
        envelopedDataDS(iVoxel,:) = decimate(HE(iVoxel,:), ...
                                             factorList(iFactor),     ...
                                             filterType);
    end%loop over parcels
    t_ds = decimate(t, factorList(iFactor), filterType);
    
    if iFactor < length(factorList), % if we're not in last iteration of loop
        HE = envelopedDataDS;
        t  = t_ds;
        clear envelopedDataDS t_ds
    end%if
end%if
end%filter_and_downsample



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [envelopedData, t_ds, newFs] = moving_average_downsample(t, overlap,   ...
                                                           windowLength, ...
                                                           useHanningWindow)
%MOVING_AVERAGE_DOWNSAMPLE downsamples signal using a moving average filter
%
% [SIGDS, TDS] = MOVING_AVERAGE_DOWNSAMPLE(SIG, T, OVERLAP, L, USEWINDOW)
%   downsamples signal SIG using a moving average filter with section length
%   L seconds, and fractional overlap OVERLAP. Each column of SIG is a
%   sample in time, rows of SIG hold different signals. T is a vector
%   indexing the time at each sample of SIG. A constant sampling rate is
%   assumed. USEWINDOW is a logical flag to apply a Hanning window to each
%   section of the moving average. 
%
%   The function returns the downsampled signal, SIGDS, and a downsampled
%   time vector, TDS. 
%
%   See also: OSL_MOVAVG. 

assert(overlap <= 1 && overlap >= 0, ...
       [mfilename ':InvalidOverlapParameter'], ...
       'The window overlap should be a fraction between 0 and 1. \n');

Fs       = 1.0 / median(diff(t));
nVoxels  = ROInets.rows(HE);
winsize  = round(windowLength * Fs);

for iVoxel = nVoxels:-1:1, % loop backwards to initialise correctly
    % moving average
    envelopedData(iVoxel,:) = osl_movavg(HE(iVoxel,:), ...
                                         t, ...
                                         winsize, ...
                                         overlap, ...
                                         useHanningWindow);
end%loop over parcels

DsFactor = length(t) / size(envelopedData, 2);
newFs    = Fs / DsFactor;
t_ds     = t(1) : 1.0/newFs : t(end);
t_ds     = t_ds(1:ROInets.cols(envelopedData));
end%moving_average_downsample





end%envelope_data






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = parse_inputs(varargin)
% parse param-value inputs to envelope_data

P               = inputParser;
P.CaseSensitive = false;
P.FunctionName  = mfilename;
P.StructExpand  = true;  % If true, can pass parameter-value pairs in a struct
P.KeepUnmatched = false; % If true, accept unexpected inputs

numericValidFcn = @(x) ( isnumeric(x) &&  isscalar(x) && ...
                        ~isnan(x)     && ~isinf(x)      );

windowLengthDefault = 2.0; % s
P.addParamValue('windowLength',     windowLengthDefault, numericValidFcn);
P.addParamValue('overlap',          0,                   numericValidFcn);
P.addParamValue('useHanningWindow', false,               @islogical);
P.addParamValue('useFilter',        false,               @islogical);
P.addParamValue('verbose',          false,               @islogical);
P.addParamValue('takeLogs',         false,               @islogical);
P.parse(varargin{:});

varargout = {P.Results.windowLength,     ...
             P.Results.overlap,          ...
             P.Results.useHanningWindow, ...
             P.Results.useFilter,        ...
             P.Results.takeLogs,         ...
             P.Results.verbose};
end%parse_inputs



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fList = get_factors(f)
% break downsample factor into smaller chunks

fMax = 13; % maximum stable downsample order

if f <= fMax,
    fList = f; 
    
else
    ff = factor(f);
    
    if     any(2 == ff) && any(5 == ff),
        fList = [10, get_factors(f ./ 10)];
        
    elseif any(2 == ff) && any(3 == ff),
        fList = [ 6, get_factors(f ./ 6)];
        
    elseif sum(2 == ff) >= 3,
        fList = [ 8, get_factors(f ./ 8)];
        
    elseif sum(2 == ff) == 2,
        fList = [ 4, get_factors(f ./ 4)];
        
    else
        fList = ff;
    end%if
end%if
end%get_factors



% [EOF]
