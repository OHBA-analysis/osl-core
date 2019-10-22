function oat = setup_source_recon_params_for_blank_recon(P)
%SETUP_SOURCE_RECON_PARAMS sets oat structure for setting up simulation
% The parameters here are optimised for running on a blank dataset
% The purpose of such a run is to extract the co-ordinate systems and lead
% fields. 


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
%   $Revision: 373 $
%   $LastChangedDate: 2015-01-12 16:36:30 +0000 (Mon, 12 Jan 2015) $
%   Contact: giles.colclough 'at' magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45


oat = struct();
% oat.paradigm = 'rest';

oat.source_recon.conditions = {'Undefined'};
oat.source_recon.time_range = [-Inf Inf];               % use all if set to ''
oat.source_recon.freq_range = [4 30];           % frequency range in Hz
[saveDir, saveFileName]     = fileparts(P.saveFileName);

oat.source_recon.dirname = fullfile(saveDir, ...
                                    [saveFileName '-oat-source-recon']); % saved location of oat structure

% oat.source_recon.sessions_to_do = [];              % set automatically
% oat.source_recon.epochinfo                    % set automatically
oat.source_recon.method       = 'beamform';
% oat.source_recon.modalities                   % set automatically from scanner type
% oat.source_recon.normalise_method             % set automatically
oat.source_recon.pca_dim      = 102;            % force pca as BF on blank data finds degenerate PCs
% oat.source_recon.regpc                        % set automatically
oat.source_recon.gridstep     = P.spatialRes;   % mm
% oat.source_recon.mask_fname                   % set automatically
% oat.source_recon.useheadshape = 1;
% oat.source_recon.use_rhino    = 0;
oat.source_recon.forward_meg  = 'MEG Local Spheres';
% oat.source_recon.beamformer_fn                % set automatically        %%% Does including this create problems in osl1.5.0? I think so - 15 May 2014
% oat.source_recon.mri          = {P.structuralFile};
% oat.source_recon.work_in_pca_subspace = 0;      % don't beamform in reduced subspace using local spheres
oat.source_recon.force_pca_dim = 1;             % force pca as BF on blank data finds degenerate PCs
% oat.source_recon.mni_coords                   % set automatically
% oat.source_recon.type                         % set automatically
% oat.source_recon.bandstop_filter_mains        % set automatically
end%setup_source_recon_params_for_blank_recon
% [EOF]