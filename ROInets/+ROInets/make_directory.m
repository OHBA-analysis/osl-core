function [] = make_directory(directoryName)
%MAKE_DIRECTORY  Makes directory if not already existing - wrapper on mkdir
%
% MAKE_DIRECTORY(NAME) checks for the existence of directory NAME and
%   creates it if it doesn't exist. 
%   Any errors are caught and displayed. 
%	
%	See also MKDIR. 

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
%	$Revision: 214 $
%	$LastChangedDate: 2014-07-24 12:40:42 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 14-Apr-2014 16:53:27

if ~isdir(directoryName), 
    [success, errMsg, errId] = mkdir(directoryName);
    if ~success, 
        error([mfilename ':' errId],               ...
              ['Creation failed of directory: \n', ...
               '   %s\n\n ',                       ...
               'Error message:\n%s\n'],            ...
              directoryName, errMsg); 
    end%if
end%if

end%make_directory