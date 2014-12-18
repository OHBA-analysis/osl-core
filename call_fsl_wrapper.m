function result = call_fsl_wrapper(fslCommand, quiet)
%CALL_FSL_WRAPPER wrapper on call_fsl function to check for errors. 
%
% RESULT = CALL_FSL_WRAPPER(COMMAND) runs an fsl function using terminal command
%    line string COMMAND, returning standard output in RESULT. 
%
% See also: CALL_FSL. (Part of the fsl distribution.)


%	Copyright 2013 OHBA
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
%	$Revision: 313 $
%	$LastChangedDate: 2014-11-11 12:29:19 +0000 (Tue, 11 Nov 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 07-Nov-2013 12:43:45
if nargin < 2 || ~exist('quiet', 'var') || ~quiet,
    fprintf(fslCommand);
    fprintf('\n');
end%if

[status, result] = call_fsl(fslCommand);
if status, 
    error([mfilename ':fslCallFailed'], ...
          'Call to fsl failed with message: \n   %s \n', result);
end