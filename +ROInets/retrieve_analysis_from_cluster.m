function [ROInetworks, correlationMats] = retrieve_analysis_from_cluster(ROInetworks,           ...
                                                                    individualSaveFileNames, ...
                                                                    deleteOriginals)
%RETRIEVE_ANALYSIS_FROM_CLUSTER(OIL) collects and formats results after
%   batching off to fmrib cluster. 
%
% RETRIEVE_ANALYSIS_FROM_CLUSTER(OIL, SAVE_FILES, DELETE_ORIGINALS) uses
%   individual session network matrices specified by full paths in cell
%   array SAVE_FILES. 
%   If DELETE_ORIGINALS is true, the session-specific matrices will be 
%   deleted from the disc.
%
% See also: osl_network_analysis, submit_analysis_to_cluster


%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 373 $
%	$LastChangedDate: 2015-01-12 16:36:30 +0000 (Mon, 12 Jan 2015) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 01-May-2014 12:17:57

fprintf('%s: Retrieving analysis results. \n', mfilename);

%% Parse input settings
fprintf('%s: checking inputs. \n', mfilename);

% oil      = ROInets.check_inputs_osl(oil); - probably unnecessary
Settings = ROInetworks;

if nargin < 2 || isempty(individualSaveFileNames),
    % use the default vales from ROInets.submit_analysis_to_cluster
    for iSession = Settings.nSessions:-1:1;
        individualSaveFileNames{iSession} = fullfile(Settings.outputDirectory,                                 ...
                                                     Settings.sessionName{iSession});
        assert(logical(exist(individualSaveFileNames{iSession}, 'file')), ...
               [mfilename ':DataNotFoundOnDisk'], ...
               'Individual network results file not found on disk: \n%s.\n', ...
               individualSaveFileNames{iSession});
    end%for
end%if

% check the files exist
if iscell(individualSaveFileNames) 
    doesExist =  cellfun(@(c) logical(exist(c, 'file')), individualSaveFileNames);
else
    error([mfilename ':NonCellFileList'], ...
          'Please input a cell array of mat files holding session network results. \n');
end%if
if ~all(doesExist),
    error([mfilename ':FilesDontExist'], ...
          'Are you sure all those files actually exist?\n');
end%if

if nargin < 3,
    deleteOriginals = false;
end%if

%% Load matrices
fprintf('%s: Loading Results. \n', mfilename);
corrMatVariableName = 'CorrMats'; % 'out' in older versions
for iSession = length(individualSaveFileNames):-1:1,
    tmp = load(individualSaveFileNames{iSession}, corrMatVariableName);
    mats{iSession} = tmp.(corrMatVariableName);
    clear tmp
end%for

% reformat results - correlationMats is a cell array of frequency bands
correlationMats = ROInets.reformat_results(mats, Settings);
clear mats

%% Group-level analysis
correlationMats = ROInets.do_group_level_statistics(correlationMats, Settings);

%% save matrices
fprintf('%s: Saving Results. \n', mfilename);
saveFileName = fullfile(Settings.outputDirectory, 'ROInetworks_correlation_mats.mat');
save(saveFileName, 'correlationMats');

if deleteOriginals,
    for iSession = 1:length(individualSaveFileNames),
        delete(individualSaveFileNames{iSession});
    end%for
end%if
 
ROInetworks.correlationMatsFile = saveFileName;

fprintf('%s: Analysis complete. \n\n', mfilename);
end%retrieve_analysis_from_cluster
% [EOF]