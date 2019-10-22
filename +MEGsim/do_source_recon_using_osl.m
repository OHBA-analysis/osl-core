function [ReconResults, oil] = do_source_recon_using_osl(oil, tmpFile, doEnv)
%DO_SOURCE_RECON_USING_OSL performs source reconstruction using OSL tools
%
% [RECONRESULTS, RECONOAT] = DO_SOURCE_RECON_USING_OSL(OAT, D) reconstructs
%   spm MEEG object D using source reconstruction settings held in
%   structure OAT. This produces the beamformer results structure
%   RECONRESULTS and the oil or oat structure, after beamforming, RECONOAT. 
% 
% [RECONRESULTS, RECONOAT] = DO_SOURCE_RECON_USING_OSL(OAT, D, DOENVELOPING) 
%   will also envelope the reconstructed data if DOENVELOPING is true. 
%   (This is the second stage of the oat pipeline.) The window length is 
%   set in oil.enveloping.window_length. 


%   Copyright 2014 OHBA
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


%   $LastChangedBy: giles.colclough@gmail.com $
%   $Revision: 213 $
%   $LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%   Contact: giles.colclough 'at' magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45

if nargin < 3 && ~exist('doEnv', 'var'), 
    doEnv = false;
end%if

% add meeg object to oil structure
D = spm_eeg_load(tmpFile);
if strcmpi(type(D), 'continuous'),
    oil.source_recon.D_continuous = {tmpFile};
    
elseif strcmpi(type(D), 'single'),
    oil.source_recon.D_epoched = {tmpFile};
    
else %unknown
    error([mfilename ':UnrecognisedFileType'], ...
          'Expecting continuous or single meeg object type. \n');
end%if

% surpress figure output
defFigVis  = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
changeBack = onCleanup(@() set(0, 'DefaultFigureVisible', defFigVis));

% check structure
% oil = osl_check_oil(oil);

% run beamformer
oil.to_do = zeros(1,6);
oil.to_do(1) = 1;

oil = osl_run_oil(oil);

% load BF results
ReconResults = load(fullfile(oil.source_recon.dirname, ...
                             oil.source_recon.results_fnames{1}));
ReconResults = ReconResults.oat_stage_results;

% do enveloping
if doEnv,
    if isfield(oil.enveloping, 'name'),
        oil.enveloping = rmfield(oil.enveloping, 'name'); 
    end%if
    
    oil.to_do = [0 1 0 0 0 0];
    oil       = osl_run_oil(oil);
end%if do enveloping

end%do_source_recon_using_osl
% [EOF]