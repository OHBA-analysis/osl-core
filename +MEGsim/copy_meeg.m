function Dout = copy_meeg(D, newname, updateHistory, removeOld)
%COPY_MEEG Copy EEG/MEG data to new files
% Allows specification of a new path, unlike SPM_EEG_COPY
%
% D = COPY_MEEG(DOLD, NEWNAME) copies DOLD to a new filename
%    NEWNAME, outputing new meeg object D. 
%
% D = COPY_MEEG(DOLD, NEWNAME, UPDATEHISTORY) is a [default: true]
%    flag to update the history of D. 
%
% D = COPY_MEEG(DOLD, NEWNAME, UPDATEHISTORY, REMOVEOLD) is a
%    [default: false] flag to delete the original object DOLD. 
%
% See also: SPM_EEG_COPY. 

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
%   $Revision:$
%   $LastChangedDate:$
%   Contact: giles.colclough 'at' magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 21-Oct-2013 
%   based on spm_eeg_copy in spm8.


if nargin<3 || ~exist('updateHistory', 'var'),
    updateHistory = true;
end%if
if nargin < 4 || ~exist('removeOld', 'var'),
    removeOld = false;
end%if

% get MEEG object
D = spm_eeg_load(D);

% get filename for the new dataset
if nargin<2 || ~exist('newname', 'var') || isempty(newname),
    newname = spm_input('New file name:', '+1', 's');
end

% remove extension
[newpath, newname] = fileparts(newname);%[spm_str_manip(S.newname, 'r') '.dat'];
newname = fullfile(newpath, newname);

% copy dataset (.mat and .dat)
Dout     = clone(D, newname);
[r, msg] = copyfile(fullfile(D.fnamedat), ...
                    fullfile(Dout.fnamedat), 'f');
if ~r
    error(msg);
end

% update history
if updateHistory
    Dout = Dout.history('spm_eeg_copy', D); 
end

save(Dout);

% cleanup
if removeOld,
    D.delete;
end%if
end%copy_meeg