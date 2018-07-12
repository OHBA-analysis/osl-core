function res = bf_inverse_mne_adaptive(BF, S)
% Computes MNE filters. 
% Wrapper for SPM around osl_inverse_mne_weights
%--------------------------------------------------------------------------

%% SPM preamble 
if nargin == 0         
    Options = cfg_entry;
    Options.tag = 'Options';
    Options.name = 'Options';        
    Options.val = {};

    Noise        = cfg_entry;
    Noise.tag    = 'Noise';
    Noise.name   = 'Measured Noise';
    Noise.val    = {};
    
    mne_adaptive      = cfg_branch;
    mne_adaptive.tag  = 'mne_adaptive';
    mne_adaptive.name = 'mne_adaptive';
    mne_adaptive.val  = {Options,Noise};      
    
    res = mne_adaptive;
    
    return
elseif nargin < 2
    error([mfilename ':ArgumentCheck'], ...
          'Two input arguments are required. \n');
end%if


%% Construct inputs
SensorData = struct('cov',      osl_cov(BF.data.D, 1), ...
                    'nSamples', BF.data.D.nsamples * BF.data.D.ntrials);
LeadFields = struct('nSources', length(S.L), ...
                    'nDims',    3,           ... % dimensions in space
                    'lf',       cat(2, S.L{:}));
   

%% Run inverse function
[W,~,L] = osl_inverse_mne_weights(SensorData, LeadFields, S.Noise, S.Options);

%% Parse outputs
res.W = W;
res.class{1}.W = W;
res.class{1}.L = L;
