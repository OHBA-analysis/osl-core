function [hmm,stats] = osl_hmm_infer(data,options)
% Runs an HMM with nstates multiple times
% [hmm,stats] = osl_hmm_infer(data,options)
%
% INPUT
% data              observations (channels x samples)
% options           structure with the training options - see documentation
%
% OUTPUT
% hmm               estimated HMM model
% stats (optional)  HMM state statistics
% AB 2014

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'HMM-MAR')))

if nargin == 1
    options = struct;
end


% Parse & remove extra options to hmmmar
if isfield(options,'Ninits')
    Ninits = options.Ninits;
    options = rmfield(options,'Ninits');
else
    Ninits = 1;
end


% Set default options
if ~isfield(options,'K')
    options.K = 8;
end
if ~isfield(options,'order')
    options.order = 0;
end
if ~isfield(options,'zeromean')
    options.zeromean = 0;
end

% Check data dimensions
if size(data,1) < size(data,2)
    data = transpose(data);
end

T = size(data,1);


% Run HMM inference with multiple initialisations
FrEn = Inf;
for i = 1:Ninits
    options.inittype='EM';
    %options.initcyc = 10;
    [hmm_new, Gamma, Xi, vpath, GammaInit, residuals, fehist] = hmmmar(data,T,options);
    % keep inference if Free Energy is lower
    if fehist(end) < FrEn
        hmm = hmm_new;
        hmm.statepath = vpath;
        FrEn = fehist(end);
    end
end 
hmm.FrEn = fehist(end);
    

% Set sampling rate
if isfield(options,'Fs')
    hmm.fsample = options.Fs;
else
    
    hmm.fsample = [];
end

    
% Output covariance matrices for MVN case
if options.order == 0
    for k = 1:hmm.K
        hmm.state(k).Cov = hmm.state(k).Omega.Gam_rate./hmm.state(k).Omega.Gam_shape;
    end
end
 
if nargout > 1
    stats = osl_hmm_stats(hmm);
end
