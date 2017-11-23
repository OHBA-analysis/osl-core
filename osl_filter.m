function D2 = osl_filter(D,fband,varargin)
	% Apply temporal filter, to MEEG or to data array
	% 
	% This function applies a temporal filter based on the contents of the two-element
	% vector 'fband'
	%
	% INPUTS
	% - D - MEEG object or n_signals x n_times matrix  - so that orientation matches D(:,:) 
	% - fband - Two-element vector of frequencies (Hz)
	% 
	% Optional
	% - fs - sampling rate, must be provided if D is a matrix
	% - order - filter order
	% - prefix - If D is an MEEG, use this prefix. Set to empty string to operate in place
	%
	% OUTPUTS
	% - If D is an MEEG 
	%		- A new MEEG object with filtered data from on the original MEEG object
	%	Otherwise
	%		- An n_signals x n_times matrix of filtered data
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

	arg = inputParser;
	arg.addParameter('order',5); 
	arg.addParameter('fs',[]); 
	arg.addParameter('prefix','f'); 
	arg.parse(varargin{:});

	if isnumeric(D)

		assert(~isempty(arg.Results.fs),'Sampling rate must be specified ')

		if fband(1) == 0
			D2 = ft_preproc_lowpassfilter(D, arg.Results.fs, fband(2), arg.Results.order, 'but','twopass','reduce');
		elseif fband(2) == inf
			D2 = ft_preproc_highpassfilter(D, arg.Results.fs, fband(1), arg.Results.order, 'but','twopass','reduce');
		elseif fband(1) < 0
			D2 = ft_preproc_bandstopfilter(D,arg.Results.fs,abs(fband),arg.Results.order, 'but','twopass','reduce');
		else
			D2 = ft_preproc_bandpassfilter(D, arg.Results.fs, fband, arg.Results.order, 'but','twopass','reduce');
		end

	else

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

		S.order = arg.Results.order;
		S.type  = 'butterworth';
		S.dir   = 'twopass';
		S.prefix = arg.Results.prefix;

		D2 = spm_eeg_filter(S);
		D2 = D2.montage('switch',montage_idx);
	
	end

