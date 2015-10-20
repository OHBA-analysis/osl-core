%         ____ _     _____    _    _   _ 
%        / ___| |   | ____|  / \  | \ | |
%       | |  _| |   |  _|   / _ \ |  \| |
%       | |_| | |___| |___ / ___ \| |\  |
%        \____|_____|_____/_/   \_|_| \_|
%
%   Group Level Exploratory Analysis of Networks                              
%
%                Adam Baker 2015
%
%
% For details on how to set up a new GLEAN analysis type "help glean_setup"
% For details on particular GLEAN settings type "help glean_check"
% To run a GLEAN analysis type "run_glean(GLEAN)"
%
% If using this toolbox please consider citing the following:
%
% "Fast transient networks in spontaneous human brain activity", Baker el al., eLife, 2014
% "Spectrally resolved fast transient brain states in electrophysiological data", Vidaurre et al., (in rev) 
% "A symmetric multivariate leakage correction for MEG connectomes", Colclough et al., NeuroImage, 2015
%
% Files:
%   glean_check               - Checks the settings and set up directory and data structures in GLEAN.
%   glean_connectivityprofile - Computes the "connectivity profile" P for each HMM state. 
%   glean_data                - Set up the directory structure for a new or existing GLEAN analysis.
%   glean_envelope            - Runs the envelope stage of GLEAN.
%   glean_hilbenv             - Optimised Hilbert envelope computation for of MEEG data. 
%   glean_model               - Runs the model stage of GLEAN.
%   glean_normalise           - Performs normalisation of MEEG data.
%   glean_parcellation        - Computes node time series from a MEEG object using a parcellation.
%   glean_regress             - Create spatial maps via mutliple regression of HMM or ICA time courses.
%   glean_run                 - Runs a GLEAN analysis.
%   glean_setup               - Sets up a GLEAN analysis with particular filename, settings and data.
%   glean_subspace            - Runs the subspace stage of GLEAN.
%   glean_convert2spm         - Converts data to SPM12 format, saving the new .mat and .dat files.
%   glean_directories         - Sets up the directory structure for the GLEAN analysis.
%   glean_groupcov            - Efficiently compute a group covariance from multiple SPM files.
%   glean_infer_hmm           - Infers an hidden Markov model (HMM) with particular options.
%   glean_plot_timecourse     - Plots GLEAN inferred time courses.
%   glean_results             - Runs the results stage of GLEAN.
