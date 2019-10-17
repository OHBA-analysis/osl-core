function MEGmodalities = find_MEG_chantypes(D)
%FIND_MEG_CHANTYPES finds MEG modalities in spm object
%
% MEGCHANTYPES = FIND_MEG_CHANTYPES(D) returns a cell array of
%   MEGCHANTYPES present in spm meeg object D. 
%   If no MEG modalities are found, an empty cell is returned. 
%

%   This code was largely borrowed from osl_check_oil.m, written by Henry
%   Luckhoo and Mark Woolrich. 

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

chanList = chantype(D, 'MEG');

if sum(ismember(unique(chanList),'MEGGRAD')),
    % ctf
    MEGmodalities = {'MEGGRAD'};
elseif sum(ismember(unique(chanList),'MEG')),
    % some other scanner type
    MEGmodalities = {'MEG'};
    
else    
    MEGmodalities = {};
    % accept neuromag sensor types -others not supported
    if ( sum(ismember(unique(chanList),'MEGPLANAR')) && ...
        ~sum(ismember(unique(chanList),'MEGPCACOMP'))),
        MEGmodalities = [MEGmodalities {'MEGPLANAR'}];
    end
    if ( sum(ismember(unique(chanList),'MEGMAG')) && ...
        ~sum(ismember(unique(chanList),'MEGPCACOMP'))),
        MEGmodalities = [MEGmodalities {'MEGMAG'}];
    end
end%if
