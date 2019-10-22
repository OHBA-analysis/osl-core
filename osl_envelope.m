function [env,phase] = osl_envelope(D,varargin)
    % Return Hilbert envelope and corresponding phases
    %
    % INPUTS
    % - D : MEEG object OR n_signals x n_samples matrix
    % - varargin : Optional inputs, see below
    %
    % OUTPUTS
    % Note that the number of output samples may be different to
    % the number of input samples if 'downsample' is specified
    %
    % - env: envelope timecourse : n_signals x n_output_samples 
    % - phase: phase timecourse : n_signals x n_output_samples 
    % 
    % NOTE - If an meeg is passed in with bad samples present, the 
    % output env and phase variables will contain NaNs at the 
    % affected times.
    %
    % Romesh Abeysuriya 2017
    
    arg = inputParser;
    arg.addParameter('fs',[]); % If D is a matrix, you must specify the sampling rate
    arg.addParameter('clean',[]); % If D is a matrix, optionally specify good samples
    arg.addParameter('filter',[]); % specify range of frequencies to filter at prior to enveloping (e.g. [8 12], default: no filter)
    arg.addParameter('downsample_env',[]); % Resample the envelope at this frequency
    arg.addParameter('downsample_phase',[]); % Resample the phase at this frequency
    arg.addParameter('orthogonalise',false); % Orthogonalise time series after filtering but before enveloping 
    arg.addParameter('orthog_method','symmetric'); % If string, can be 'symmetric' or 'closest'. Or, use a number to specify PM filter order
    arg.parse(varargin{:});

    if isa(D,'meeg')
        fs = D.fsample;
        clean = good_samples(D);
        D = D(:,:);
    else
        fs = arg.Results.fs;
        if isempty(fs)
            error('If D is a matrix, the sampling rate must be specified to ensure filtering and downsample work properly');
        end
        if isempty(arg.Results.clean)
            clean = ones(1,size(D,2));
        else
            clean = arg.Results.clean;
        end
    end

    % If doing PM orthogonalization, do it before filtering
    if arg.Results.orthogonalise && isnumeric(arg.Results.orthog_method)
        D = leakcorr(D.',size(D,2),arg.Results.orthog_method).';
    end

    if ~isempty(arg.Results.filter)
        D = osl_filter(D,arg.Results.filter,'fs',fs);
    end

    % If doing Giles orthogonalization, do it after
    if arg.Results.orthogonalise && ~isnumeric(arg.Results.orthog_method)
        D = ROInets.remove_source_leakage(D,arg.Results.orthog_method);
    end

    % Compute the envelope and phase timecourses
    analytic = hilbert(bsxfun( @minus, D.', mean(D.') ));
    env = abs(analytic).';
    phase = unwrap(angle(analytic)).';

    
    if ~isempty(arg.Results.downsample_env)
        env = ds_clean(env,clean,fs,arg.Results.downsample_env);
    else
        env(:,~clean) = NaN;
    end

    if ~isempty(arg.Results.downsample_phase)
        phase = ds_clean(phase,clean,fs,arg.Results.downsample_phase);
    else
        phase(:,~clean) = NaN;
    end


function D2 = ds_clean(D,clean,f_original,f_new)
    % downsample and clean
    % Inputs
    % - D original data
    % - clean - logical vector specifying whether time is clean or not
    % - f_original - original sampling rate
    % - f_new - new sampling rate
    % Output will be resampled, and will contain NaNs for times where 
    % the original sample maps to the downsampled data (see map_times() below)

    assert(size(D,2)==length(clean));
    t_original = (0:size(D,2)-1)/f_original;
    D2 = resample(D.',f_new,f_original).';
    t_new = (0:size(D2,2)-1)/f_new;
    mapping = map_times(t_original,t_new);
    D2(:,mapping(~clean)) = NaN;

function mapping = map_times(t1,t2)
    % Given a time vector like
    % t1 = [1 2 3 4 5 6 7 8 9 10]
    % And downsampled times like
    % t2 = [1 4 8]
    % Find which index in t2 is closest to t1
    % e.g
    % mapping = [1 1 2 2 2 3 3 3 3 3]

    assert(isvector(t1) && isvector(t2),'Inputs must be vectors');
    assert(issorted(t1) && issorted(t2),'Inputs must be sorted');

    ptr = 1;
    mapping = zeros(size(t1));

    for j = 1:length(t1)
        d = abs(t1(j) - t2([ptr ptr+1]));
        if d(1) >= d(2)
            ptr = ptr + 1;
        end

        if ptr == length(t2)
            mapping(j:end) = ptr;
            break
        else
            mapping(j) = ptr;
        end

    end
