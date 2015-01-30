function D = oslview_filter(D)

fspec = inputdlg({'F1','F2'},'Filter',1,{'1','48'});

if isempty(fspec)
  return
end

fbands = [str2double(fspec{1}) str2double(fspec{2})];
fbands = sort(fbands);

S=[];
S.D    =  fullfile(D.path,D.fname);
S.freq = fbands;
S.band = 'bandpass';
Dnew = spm_eeg_filter(S);

D = clone(Dnew,fullfile(D.path,D.fname));
copyfile(Dnew.fnamedat,D.fnamedat, 'f');

save(D);
delete(Dnew);


end