function megInd = megchannels(D, lookFor)
%MEGCHANNELS returns indices of meeg channels
%
% I = MEGCHANNELS(D) returns scalar indices I of meeg channels in MEEG
%   object D. 
%
%   This replaces MEEGCHANNELS(D, 'MEG'), which misses out MEGPLANAR
%   options. 
%
% See also: MEEGCHANNELS. 


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
%   $Revision: 293 $
%   $LastChangedDate: 2014-11-05 13:32:02 +0000 (Wed, 05 Nov 2014) $
%   Contact: giles.colclough 'at' magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 29-Nov-2013 17:38:34

% default behaviour: look for meg channels
if 1 == nargin, 
    lookFor = 'MEG';
end

megInd = find(strncmpi(lookFor, chantype(D), length(lookFor)));