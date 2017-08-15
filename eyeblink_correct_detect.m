function blinkfname=eyeblink_correct(S)

disp('Loading data.')

fname = S.D;

Dorig = spm_eeg_load(fname);

%1.detect eyeblinks
disp('Detecting eyeblinks.');
S = [];
S.D = Dorig;
S.method = 'events';
S.doplots = 1;
S.sd_thr = 3;
S.fname=fname;
D = detect_eyeblinks(S);
%close all;

%2.epoch w.r.t. eyeblink
disp('Epoching w.r.t. eyeblinks.');
S = [];
S.D = D;
S.inputformat = [];
S.pretrig = -300;
S.posttrig = 300;
S.trialdef(1).conditionlabel = 'blink';
S.trialdef(1).eventtype = 'eyeblink';
S.trialdef(1).eventvalue = 1;
S.reviewtrials = 0;
S.save = 0;
S.epochinfo.padding = 0;
S.bc = 2;
S.prefix = 'blink_e'; % note this requires a modified version of spm_eeg_epochs that accepts this input argument... LH 280509

D = spm_eeg_epochs_v2(S);
blinkfname=fullfile(D.path,D.fname);

S = [];
S.D = D;

%D = spm_eeg_average(S);

%3.perform SVD
S = [];
S.D = D;
S.method = 'SVD';
S.timewin = [-300 300]/1000; %600 ms timewindow;
S.ncomp = 4;

[D,svdoutput] = spm_eeg_spatial_confounds_v2(S);

D.save;

mkdir([fname '_svdoutputs']);
svdfname = fullfile([fname '_svdoutputs'],'blink_svd.mat'];
save(svdfname,'svdoutput');




    
    
    