function newNames = prefix(oldNames, prefix)
%PREFIX  prefixes a character onto a set of filenames. 
% NEW = PREFIX(OLD, PRE) adds prefix PRE to filename OLD to make
%   NEW filename. OLD can be a single string, or a cell array of strings.
%   PRE is a character string. OLD can contain directories, these are taken
%   care of.
%   

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


%   $LastChangedBy$
%   $Revision$
%   $LastChangedDate$
%   Contact: giles.colclough@magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 24-Oct-2014 10:44:58

% input checking
assert(ischar(prefix),                     ...
       [mfilename ':IncorrectPrefixType'], ...
       'Prefix should be a string. \n');
   
% apply prefix, catering for character or cell arrays
if iscell(oldNames), 
    newNames = cellfun(@(name) add_prefix(name, prefix), ...
                       oldNames,                         ...
                       'UniformOutput', false);
    
elseif ischar(oldNames),
    newNames = add_prefix(oldNames, prefix);
    
else 
    error([mfilename ':InconsistentInput'], ...
          'Expecting a string or cell array of strings as input. \n');
end%if
end%prefix

function new = add_prefix(old, pre)
%ADD_PREFIX

[dirName, fName, ext] = fileparts(old);

new = fullfile(dirName, sprintf('%s%s%s', pre, fName, ext));
end%add_prefix
% [EOF]