function D2 = osl_filter(D,fband,order)
	% Apply temporal filter, return new MEEG object
	% 
	% This function applies a temporal filter based on the contents of the two-element
	% vector 'fband'
	%
	% INPUTS
	% - D - MEEG object
	% - fband - Two-element vector of frequencies (Hz)
	% - order (optional, otherwise will use spm_eeg_filter's default)
	%
	% OUTPUTS
	% - A new MEEG object with filtered data from on the original MEEG object
	%
	% fband both selects the filter type and the frequencies. Suppose the input is
	% 
	% fband = [f1 f2]
	%
	% Options are:
	% 
	% f1 = 0 - Low pass filter with cutoff f2
	% f2 = inf - High pass filter with cutoff f1
	% f1 < 0 - Bandstop filter with range abs(f1) to abs(f2)
	% Otherwise - Bandpass filter with range f1 to f2
	% 
	% EXAMPLES
	% 
	% D2 = osl_filter(D,[0 10]) % Low-pass filter with cutoff at 10Hz
	% D2 = osl_filter(D,[45 inf]) % High-pass filter with cutoff at 45Hz
	% D2 = osl_filter(D,-[48 52]) % Bandstop filter blocking frequencies 48-52Hz
	% D2 = osl_filter(D,[8 13]) % Bandpass filter allowing 8-13Hz

	if nargin < 3 || isempty(order) 
		order = 5;
	end
	
	S = struct;

	montage_idx = D.montage('getindex');
	S.D = D.montage('switch',0);

	if fband(1) == 0
		S.band = 'low';
		S.freq = fband(2);
	elseif fband(2) == inf
		S.band = 'high';
		S.freq = fband(1);
	elseif fband(1) < 0
		S.band = 'stop';
		S.freq = abs(fband);
	else
		S.band = 'bandpass';
		S.freq = fband;
	end

	S.order = order;
	S.type  = 'butterworth';
	S.dir   = 'twopass'

	D2 = spm_eeg_filter(S);
	D2 = D2.montage('switch',montage_idx);
	



