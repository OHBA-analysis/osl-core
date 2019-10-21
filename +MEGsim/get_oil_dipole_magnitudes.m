function OILdipoleMag = get_oil_dipole_magnitudes(simDataReconResults, ...
                                                  iTrial)
%GET_OIL_DIPOLE_MAGNITUDES extract dipole magnitudes from source recon


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

oilResults    = MEGsim.AB_get_source_timecourses(simDataReconResults);

% just use first trial
turn_into_mag = @(C) sum(C(:,:,iTrial).^2, 1);
OILdipoleMagCell = cellfun(turn_into_mag, ...
                           oilResults.source_timecourses, ...
                           'UniformOutput', false);
% we assume HMM stuff not being done. This is true for the settings
% in this file. Otherwise the indexing below may well break.
assert(size(OILdipoleMagCell, 2)==1);

OILdipoleMag = cat(1, OILdipoleMagCell{:});
end%get_oil_dipole_magnitudes