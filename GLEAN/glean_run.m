function glean_run(GLEAN)                                 
% Top level function for running GLEAN 
%
%         ____ _     _____    _    _   _ 
%        / ___| |   | ____|  / \  | \ | |
%       | |  _| |   |  _|   / _ \ |  \| |
%       | |_| | |___| |___ / ___ \| |\  |
%        \____|_____|_____/_/   \_|_| \_|
%
%   Group Level Exploratory Analysis of Networks                              
%
%                  GLEAN v0.1 
%                Adam Baker 2015
%
% OSL GLEAN
% OSL Group Level Exploratory Analysis of Networks

% FILE STRUCTURE:
% 
% - ENVELOPES_{settings}
%   - data
%   - SUBSPACE_{settings}
%     - data
%     - model_{settings}
%       - RESULTS_{settings}
%         - PLOTS_{settings}
%
% GLEAN STRUCTURE:
%
% SETTINGS
% DATA
% MODEL
% OUTPUT
%

% TODO
% - allow running the GLEAN using the GLEAN
% - sort out input parsing & logic
% - add multiband support
% - add ability to switch montages
% - add more results options (conn. profile, stats, graphs etc)
% - normalisation of PCs/parcels pre-concatentation
% - write some default settings for the usual pipelines
%   (eLife,parcellation,ica)
% - If new envelope files are added then it needs to redo everything!
% - Add some helpful command line print outs
% - Consider moving the filename set up to the individual stages
% - Make single frequency mode a special case of multiband (i.e. make TF)

% Either create a new GLEAN or load an existing one and modify its
% parameters

save(GLEAN.name,'GLEAN')   

% COMPUTE ENVELOPES
glean_envelope(GLEAN)

% COMPUTE SUBSPACE DATA
glean_subspace(GLEAN)

% RUN THE NETWORK MODEL (ICA/HMM)
glean_model(GLEAN)

% CREATE OUTPUT MAPS/STATS
glean_results(GLEAN)



end






                    
