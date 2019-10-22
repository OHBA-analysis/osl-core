function [] = check_filename_is_clear(filename)
%CHECK_FILENAME_IS_CLEAR checks for and removes existing objects
%
% CHECK_FILENAME_IS_CLEAR(FILENAME) checks for a MEEG object at
%   FILENAME. If one exists, it will be overwritten later - delete
%   the existing object now.
%   If a file is there, but it is not a MEEG object, error in case
%   of accidental over-writes.
%   If nothing is there, proceed.


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

if osl_util.isfile(filename) || osl_util.isfile([filename '.mat']),
    try
        % delete existing object, so we can replace cleanly later
        Dtest = spm_eeg_load(filename);
        
        warning([mfilename ':OverwritingPreviousObject'], ...
                'The MEEG object at %s will be overwritten. \n', ...
                filename);
            
        successfulDeletion = delete(Dtest);
%         if ~successfulDeletion,
%             error([mfilename ':FailedToDeleteExistingObject'], ...
%                   'Failed to delete existing object %s. \n', ...
%                   filename);
%         end%if
%         clear Dtest
    catch ME 
        % check that filename wasn't an spm object
        if strcmp(ME.message, ...
                 sprintf('Trouble reading file %s', filename)),
             error([mfilename ':OverwritingUnknownFile'], ...
                   ['The filename provided referred to a non-MEEG object, \n', ...
                   '%s. \n  To prevent files being over-written, please ', ...
                   'use another file name. \n'], ...
                   filename);
        else % it was an spm object, and a different error occured
            rethrow(ME);
        end%if not spm object 
    end%try
end%if filename already exists
end%check_filename_is_clear
