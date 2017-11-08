%% HOW TO MAXFILTER

fif_file = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif';

%% First, do a nosss
nosss_fif = osl_maxfilter(fif_file,'nosss');

%% Then, find bad channels
D = osl_import(nosss_fif);
D = osl_detect_artefacts(D,'badtimes',false);

%% Now use these to maxfilter
% Don't forget to also capture bad times...
[sss_fif,bad_segments] = osl_maxfilter(fif_file,'sss','badchannels',D.chanlabels(D.badchannels));

%...and include them when importing
D_sss = osl_import(sss_fif,'bad_segments',bad_segments);

% Did that fail? How aout tsss?
[sss_fif,bad_segments] = osl_maxfilter(fif_file,'sss','badchannels',D.chanlabels(D.badchannels));
D_tsss = osl_import(sss_fif,'bad_segments',bad_segments);


%%
%% raw data
S = [];
S.fif_file = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif';
S.spm_file = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard_raw';
D = osl_import( S );

%% nosss

S = [];
S.fif = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif';
S.logfile = 1;
S.fif_out = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/nosss_fif_spm_meg1';
S.nosss = 1;
S.downsample_factor = 4; % optional
fif_out = osl_call_maxfilter(S);

S2 = [];
S2.fif_file = [S.fif_out '.fif'];
S2.spm_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/nosss_fif_spm_meg1';
D2 = osl_import( S2 );

D3 = osl_detect_artefacts(D2,'badtimes',false);

%% sss - no bad channels

S = [];
S.fif = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif';
S.logfile = 1;
S.fif_out = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1_nobadchannels';
S.movement_compensation = 1;
S.autobad_off = 1; %default

fif_out = osl_call_maxfilter(S);

%%
S = [];
S.fif_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1_nobadchannels.fif';
S.spm_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1_nobadchannels';
D4 = osl_import(S);

%% sss - bad channels

S = [];
S.fif = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif';
S.logfile = 1;
S.fif_out = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1';
S.spmfile= D3;
S.movement_compensation = 1;
S.autobad_off = 1; %default

fif_out = osl_call_maxfilter(S);

hpi = osl_headpos('/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1_headpos.txt',true);

%%
S = [];
S.fif_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1.fif';
S.spm_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1';
D5 = osl_import(S);


%% tsss

S = [];
S.fif = '/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif';
S.logfile = 1;
S.fif_out = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/tsss_fif_spm_meg1';
S.spmfile= D3;
S.movement_compensation = 1;
S.autobad_off = 1; %default
S.temporal_extension = 1;

fif_out = osl_call_maxfilter(S);

S = [];
S.fif_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/tsss_fif_spm_meg1.fif';
S.spm_file = '/home/disk3/ajquinn/Projects/lizard/first_pass.opt/tsss_fif_spm_meg1';
D5 = osl_import(S);

%%

'mne_show_fiff --in /home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1.fif --blocks | grep SSS'

hdr = ft_read_header('/home/disk3/ajquinn/Projects/lizard/raw_data/s01_lizard.fif','checkmaxfilter',true);
is_maxfilter = isempty(hdr.orig.raw.info.projs)

hdr = ft_read_header('/home/disk3/ajquinn/Projects/lizard/first_pass.opt/sss_fif_spm_meg1.fif','checkmaxfilter',true);
is_maxfilter = isempty(hdr2.orig.raw.info.projs)
