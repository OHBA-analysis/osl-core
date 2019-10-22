function y = osl_surrogate(x,method,varargin)
    % Thin wrapper to generate surrogate data using various methods
    %
    % INPUTS
    % - x : Matrix of data, nsamples x nchannels
    % - method : Method used for generate surrogate data. Options include
    %   - 'phase_randomization' : Fourier phase randomization (default)
    %   - 'ar' : mass univariate autoregressive model
    %   - 'mar' : multivariate autoregressive model
    % - varargin : passed to randomization function e.g. 'ar' and 'mar' require
    % that the model order is specified
    %
    % OUTPUTS
    % - y : Surrogate data corresponding to real input data 'x'
    % 
    % EXAMPLE USAGE
    %   y = osl_surrogate(x,'phase_randomization');
    %   y = osl_surrogate(x,'phase_randomization',true,true);
    %   y = osl_surrogate(x,'ar',10);


    if nargin < 2 || isempty(method) 
        method = 'phase_randomization';
    end
    
    assert(ndims(x)==2,'Incorrect data dimension, must be a matrix')
    if size(x,2)>size(x,1)
        fprintf(2,'Warning - data has more columns than rows. osl_surrogate expects rows to correspond to samples\n')
    end

    switch method
        case 'phase_randomization'
            y = ROInets.generate_phase_surrogates(x, varargin{:});
        case 'ar'
            y = generate_ar_surrogates(x, varargin{:});
        case 'mar'
            % Third input to osl_surrogate() is expected to be the model order
            hmm = hmmmar(x,size(x,1),struct('K',1,'order',varargin{1},'updateGamma',0,'covtype','full','verbose',0));
            y = simhmmmar(size(x,1),hmm);
        otherwise
            error('Unrecognized surrogate method');
    end

end
