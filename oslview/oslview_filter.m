function D = oslview_filter(D)

fspec = inputdlg({'F1','F2'},'Filter',1,{'1','48'});

if isempty(fspec)
  return
end

fbands = [str2double(fspec{1}) str2double(fspec{2})];
fbands = sort(fbands);

S=[];
S.D=D;
S.filter.PHz=fbands;
S.filter.band='bandpass';
Dnew = spm_eeg_filter_v2(S);


D = clone(Dnew,fullfile(D.path,D.fname));
copyfile(fullfile(Dnew.path, Dnew.fnamedat),fullfile(D.path, D.fnamedat), 'f');

save(D);
delete(Dnew);


end