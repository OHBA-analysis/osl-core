function [env,phase] = osl_envelope(D,varargin)
	% Return Hilbert envelope and corresponding phases
    %
    % INPUTS
    % - D : MEEG object OR n_samples x n_signals matrix
    % - varargin : Optional inputs, see below
    %
    % OUTPUTS
    % Note that the number of output samples may be different to
    % the number of input samples if 'downsample' is specified
    %
    % - env: envelope timecourse : n_output_samples x n_signals
    % - phase: phase timecourse : n_output_samples x n_signals
    %
    % Romesh Abeysuriya 2017
    
    arg = inputParser;
    arg.addParameter('fs',[]); % If D is a matrix, you must specify the sampling rate
    arg.addParameter('filter',[]); % specify range of frequencies to filter at prior to enveloping (e.g. [8 12], default: no filter)
    arg.addParameter('downsample_env',[]); % Resample the envelope at this frequency
    arg.addParameter('downsample_phase',[]); % Resample the phase at this frequency
    arg.addParameter('orthogonalize',false); % Orthogonalize time series after filtering but before enveloping 
    arg.addParameter('orthog_method','symmetric'); % If string, can be 'symmetric' or 'closest'. Or, use a number to specify PM filter order
    arg.parse(varargin{:});

    if isa(D,'meeg')
    	fs = D.fsample;
    	D = D(:,:);
    else
    	fs = arg.Results.fs;
    end

    % If doing PM orthogonlization, do it before filtering
    if arg.Results.orthogonalize && isnumeric(arg.Results.orthog_method)
    	D = leakcorr(D.',size(D,2),arg.Results.orthog_method).';
	end

    if ~isempty(arg.Results.filter)
        D = osl_filter(D,arg.Results.filter,[],fs);
    end

    % If doing Giles orthogonalization, do it after
    if arg.Results.orthogonalize && ~isnumeric(arg.Results.orthog_method)
    	D = ROInets.remove_source_leakage(D,arg.Results.orthog_method);
    end
	  
    % Compute the envelope and phase timecourses
    analytic = hilbert(bsxfun( @minus, D.', mean(D.') ));
    env = abs(analytic).';
    phase = unwrap(angle(analytic)).';

    if ~isempty(arg.Results.downsample_env)
    	env = resample(env,arg.Results.downsample_env,fs)
    end

    if ~isempty(arg.Results.downsample_phase)
    	phase = resample(phase,arg.Results.downsample_phase,fs)
    end


