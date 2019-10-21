function oil = setup_source_recon_params(P)
%SETUP_SOURCE_RECON_PARAMS sets oil structure for source recon pipeline
% The parameters here are optimised for running on a simulated data set,
% using nearly identical parameters to those used to create the simulation.
% 
% The only difference is that the pca dimensionalities can be left free. 


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

% use mostly the same settings from the blank recon
% -- the whole point is to beamform using the same settings as for
%    simulation
oil = MEGsim.setup_source_recon_params_for_blank_recon(P);

% We will want to change the pca properties
oil.source_recon.pca_dim       = 0;
oil.source_recon.force_pca_dim = false; % don't need to force this to happen

% Set enveloping properties
oil.enveloping.window_length = 2; %s

end%setup_source_recon_params
% [EOF]