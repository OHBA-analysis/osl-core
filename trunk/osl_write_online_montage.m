function D = osl_write_online_montage(S)
% Writes the current online montage to a new MEEG object (with prefix 'M')
% FORMAT D = osl_write_online_montage(S)
%
% S                     - input structure 
%  fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%                         continuous data
%
% Output:
% D                     - MEEG object (also written on disk)
%
% Adam Baker 2015

D = spm_eeg_load(S.D);

S               = [];
S.D             = fullfile(D.path,D.fname);
S.mode          = 'write';
S.montage       = montage(D,'getmontage');
S.keepothers    = 0;
S.keepsensors   = 0;
S.updatehistory = 1;
S.prefix        = 'M';

D = spm_eeg_montage(S);

end