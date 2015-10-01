function GLEAN = glean_setup(glean_name,data,settings)
% Sets up a GLEAN analysis with particular filename, settings and data
%
% glean_name:
%
% - full path and filename for the GLEAN analysis. The metadata and results 
%   will be saved in a .mat file by this name. Any data created by GLEAN 
%   will be saved in subfolders within the same directory as filename. 
%   Note that it is possible to create multiple GLEANs with different names 
%   within the same directory. GLEAN analyses may be run with different 
%   settings may be contained in the same directory but the data must NOT 
%   change.
%
% data:
%
% - list of OSL beamformed data files
%
% settings:
%
% - list of settings for the GLEAN analysis with the following subfields:



% GLEAN name
[pathstr,filestr,extstr] = fileparts(glean_name);
if isempty(extstr)
    extstr = '.mat';
elseif ~strcmp(extstr,'.mat')
    error('glean_name extension should be .mat');
end

if isempty(pathstr)
    error('glean_name should specify the full path and filename')
end

GLEAN.name = fullfile(pathstr,[filestr,extstr]);
GLEAN.data = data;

% Settings:
%GLEAN.settings = glean_settings(settings);
GLEAN.envelope.settings = settings.envelope;
GLEAN.subspace.settings = settings.subspace;
GLEAN.model.settings    = settings.model;
GLEAN.results.settings   = settings.results;

% Data:
GLEAN = glean_data(GLEAN);




















