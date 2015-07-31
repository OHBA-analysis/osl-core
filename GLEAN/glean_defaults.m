function GLEAN = glean_defaults(method)
% !!! THIS FUNCTION IS UNDER CONSTRUCTION !!!
% 
% Sets default settings for GLEAN if using ICA or HMM. 
%
% GLEAN = glean_defaults(method)
% 
% methods may be 'ica' or 'hmm'
% 
% 
% 
% switch method
%     
%     case 'ica'
%         
%         % Envelope settings
%         envelope.log                            = 0;
%         envelope.winsize                        = 2;
%         envelope.normalisation                  = 'global';
%         envelope.weights_normalisation          = 1;
%         
%         % Subspace settings
%         subspace.pca.dimensionality             = 20;
%         subspace.pca.whiten                     = 1;
% 
%         % Decompositon settings
%         model.ica.order = 20;
%         
%         % Output settings
%         output.pcorr.format = 'nii';
%         
%     case 'hmm'
%         
%         % Envelope settings
%         envelope.log                            = 0;
%         envelope.winsize                        = 0.1;
%         envelope.normalisation                  = 'global';
%         envelope.weights_normalisation          = 1;
%         
%         % Subspace settings
%         subspace.pca.dimensionality             = 40;
%         subspace.pca.whiten                     = 1;
%         
%         % Decompositon settings
%         model.hmm.nstates = 8;
%         model.hmm.nreps   = 3;
%         
%         % Output settings
%         output.pcorr.format = 'nii';
% 
% end