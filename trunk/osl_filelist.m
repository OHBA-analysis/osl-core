function sortedList = osl_filelist(pathstr,filestr)
%OSL_FILELIST Gets a list of files from a directory matching a pattern
%
% LIST = OSL_FILELIST(DIRNAME, PATTERN) matches files with string PATTERN
%   in directory DIRNAME, sorting in a sensible order. 
%
%
% Adam Baker 2014


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

%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$


tmp        = dir(fullfile(pathstr,filestr));
filelist   = {tmp(:).name};
filelist   = fullfile(pathstr,filelist(:));
sortedList = sort_nat(filelist);
    
end%osl_filelist