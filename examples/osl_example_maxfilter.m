%% HOW TO MAXFILTER

fif_file = 's01_lizard.fif';

%% Basic usage

fif_out = osl_maxfilter(fif_file,'nosss'); % Read only
fif_out = osl_maxfilter(fif_file,'sss'); % Do SSS
fif_out = osl_maxfilter(fif_file,'tsss'); % Do temporal SSS

%% Double maxfilter
% Standard procedure

%% 
% First, do a nosss
nosss_fif = osl_maxfilter(fif_file,'nosss');

% Then, find bad channels
D = osl_import(nosss_fif);
D = osl_detect_artefacts(D,'badtimes',false);

%%
% Now use these to maxfilter
[sss_fif,bad_segments] = osl_maxfilter(fif_file,'sss','badchannels',D.chanlabels(D.badchannels),'fif_out','awesome_sss_final.fif','verbose',true);

%%
% Note that bad times get set to zero, and these times are checked against the output log file. These
% times can be specified in |osl_import| so they are correctly treated as artefacts
D_sss = osl_import(sss_fif,'bad_segments',bad_segments);

%%
% Could do exactly the same thing with tsss
[sss_fif,bad_segments] = osl_maxfilter(fif_file,'tsss','badchannels',D.chanlabels(D.badchannels));
D_tsss = osl_import(sss_fif,'bad_segments',bad_segments);

