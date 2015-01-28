function D = osl_write_online_montage(D)
% Writes the current online montage to a new MEEG object (with prefix 'M')

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