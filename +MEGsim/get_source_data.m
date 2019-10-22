function [sourceData, time, Fs] = get_source_data(oat, resultsFileName, mask)
%GET_SOURCE_DATA  retrieves source-space time-courses
%
% [SOURCEDATA, TIME, FS] = GET_SOURCE_DATA(OAT, RECONRESULTSFILENAME)
%   retrieves voxels x time SOURCEDATA, together with a vector of times
%   TIME and the sampling rate FS. Bad channels and bad time segments are
%   ignored or removed. 
%
%   Data are retrieved from an OAT analysis structure, for the session with
%   recon results saved in RECONRESULTSFILENAME. This will be one of the
%   entries stored in OAT.source_recon.results_fnames.
%
% [] = GET_SOURCE_DATA(OAT, RECONRESULTSFILENAME, MASK) uses binary MASK to
%   select a subset of the data. MASK should have sufficient elements to
%   match the resolution of the source reconstruction. 
%
% [] = GET_SOURCE_DATA(RECONRESULTS, 1, MASK) uses pre-loaded OAT recon
%   results structure RECONRESULTS. 
%
%   See also: get_voxel_recon_timecourse,
%             get_voxel_recon_timecourse_vector, 
%             MEGsim.AB_get_source_timecourses. 


%   Copyright 2014 Giles Colclough
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


%   $LastChangedBy$
%   $Revision$
%   $LastChangedDate$
%   Contact: giles.colclough@eng.ox.ac.uk
%   Originally written on: MACI64 by Giles Colclough, 15-Oct-2014 12:20:55

fprintf('Extracting timecourses. \n');

% parse mask input
if nargin > 2 && exist('mask', 'var') && ~isempty(mask),
    validateattributes(mask,                   ...
                       {'logical', 'numeric'}, ...
                       {'binary', 'vector'},   ...
                       mfilename,              ...
                       'mask',                 ...
                       3);
    voxelInd = find(mask);
    argList = {'voxelind',voxelInd};
else
    argList = {};
end%if

if islogical(resultsFileName) || isnumeric(resultsFileName),
    ReconResults = oat;
    clear oat
else
    ReconResults = osl_load_oat_results(oat, resultsFileName);
end%if
SourceSpace  = MEGsim.AB_get_source_timecourses(ReconResults, argList{:});
clear ReconResults

% ABAS: since AB_get_source_timecourses rejects bad epoch unless flagged 
% otherwise, time should be defined from ABReconResults! 
time = SourceSpace.time(SourceSpace.time_inds); 
Fs   = SourceSpace.fsample;



sourceData = cat(1, SourceSpace.source_timecourses{:});
end%get_source_data
% [EOF]