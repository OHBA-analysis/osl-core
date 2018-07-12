% In this practical we will work with a single subject's EEG data
% 
% 

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

%OSLDIR = getenv('OSLDIR');
    
%tilde='/home/mwoolrich/Desktop';
tilde='/Users/robert/';

osldir=[tilde '/Applications/osl/osl2'];    

addpath(osldir);

osl_startup();

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

testdir=[tilde '/Documents/workshop/example/faces_subject1_data'];

workingdir=[testdir '/preprocessing']; % this is the directory the SPM files will be stored in

cmd = ['mkdir ' workingdir]; system(cmd); % make dir to put the results in
 
% Set up the list of subjects and their structural scans for the analysi 
% Currently there is only 1 subject.
clear fif_files spm_files structural_files

% list of fif files (we only have one here - note that it has already been
% Maxfiltered and downsampled to 250Hz)
fif_files{1}=[testdir '/fifs/sub1_face_sss.fif']; 

% set up a list of SPM MEEG object file names (we only have one here)
spm_files{1}=[workingdir '/spm8_eeg1.mat'];
spm_files_filtered{1}=[workingdir '/spm8_eeg1.mat'];
spm_files_epoched{1}=[workingdir '/espm8_eeg1.mat'];

% structural files in the same order as spm_files or fif_files:
structural_files{1}=[testdir '/meeg_data/mr_images/struct_S01.nii'];

trig_files{1}=['/Users/woolrich/vols_data/meeg_data/beh_data/S01_EEG/onerhy_pilot1_EEG'];

cleanup_files=0; % flag to indicate that you want to clean up files that are no longer needed

%%%%%%%%%%%%%%%%%%%%
%% CONVERT FROM FIF TO AN SPM MEEG OBJECT:
% The fif file that we are working with is sub1_face_sss.fif. This has
% already been maxfiltered for you and downsampled to 250Hz.
% 
% This will produce a histogram plot showing the number of events detected
% for each code on the trigger channel. The codes used on the trigger
% channel for this experiment were: 
%

for i=1:length(fif_files), % iterates over subjects
    
    S2.fif_file=fif_files{i};
    S2.spm_file=spm_files{i};
    S2.trigger_channel_mask='0000000000111111';
    
    % The conversion to SPM will show a histogram of the event codes
    % and correspond to those listed below in the epoching section
    D = osl_import(S2);
    close;
    Fs=fsample(D);
    
    % load trigger information and save it in the D structure
    
    nblocks = 6; % behavoural blocks per session 
    ntriggers=170;
    evtimes = nan(ntriggers,nblocks);
    evlabels = nan(ntriggers,nblocks);
    evvals = nan(ntriggers,nblocks);
    
    for ipart = 1:nblocks
        lab=['0' num2str(ipart)];
        lab=lab(end-1:end);
        load([trig_files{i} '_block' lab '_cfg.mat']);
        evtimes(:,length(1:ipart)) = cfg.trl(1)+cfg.evt(:,1);
    end
    evtimes = evtimes(:);
    
    for ipart = 1:nblocks;
        evlabels(:,ipart) = (sequence(ipart).tilt(:)+8) ;
        for tr = 1:length(evlabels)
            if isinf(evlabels(tr,ipart)) == 1;
                evlabels(tr,ipart)=17;
            end
        end
        evlabels(:,ipart) = evlabels(:,ipart) + (sequence(ipart).njitter(:)*10);
    end
    
    evlabels = evlabels(:);
    
    
    ev_org=events(D);
    D=events(D,[],[]);  % delete all events
    
    % LOOP TO SAVE EVENTS FROM EPOCH DATA INTO THE CONTINUOUS ONE!!!
    ev_new=[];
    for c = 1:length(evtimes)
        ev_new(c).type = 'STI101_up';
        ev_new(c).value = num2str(evlabels(c));
        ev_new(c).duration = ev_org(1).duration;
        ev_new(c).time = evtimes(c)/1000;
        ev_new(c).offset = ev_org(1).offset;
        
    end
    D=events(D,[],ev_new);
    D.save;
   
    % look at the data view OSLVIEW

    oslview_eeg(D);
    
    % remove MEG signals from the file  (just to save disk space and speed up
    % computations)
    
    
    labelorg = D.chanlabels;
    del_list=find(strncmp(chantype(D),'MEG',3));
    
    labelnew =labelorg;
    labelnew(del_list)=[];
    
    tra = eye(D.nchannels);
    tra(del_list, :) = [];
    
    S = [];
    S.D = D;
    S.montage.labelorg=labelorg;
    S.montage.labelnew=labelnew;
    S.montage.tra=tra;
    S.keepothers = 'no';
    D2 = spm_eeg_montage(S);
    
    
    % detect and remove bad channels
    
    S = [];
    S.D = D2;
    S.modality='EEG';
    S.rel_thres=4;
    [list_bch]=osl_detect_badchannel(S);
    bad_list=badchannels(D2);
    vect_bch=zeros(1,size(D2,1));
    vect_bch(bad_list)=1;
    vect_bch(list_bch)=1;
    D2=badchannels(D2,[],vect_bch);
    D2.save;
    
    
    % band-pass filter signals
    
    band_pass=[1 80];
    S = [];
    S.D = D2;
    S.dir = 'twopass';
    S.band = 'bandpass'; % !! onepass as this is only for visualization
    S.freq  = band_pass;
    S.use_fft_bandpass = 0;
    D3 = spm_eeg_filter(S);
    
    oslview_eeg(D3);
    
    
    Fs_new=200;
    S = [];
    S.D = D3;
    S.fsample_new = Fs_new;
    D4 = spm_eeg_downsample(S);
    
    
    %%
    D4=spm_eeg_load([workingdir '/dfMspm8_eeg1']);
    
    % artifact attenuation using AFRICA
    
    fname=fnamedat(D4);
    S=[];
    S.fname=[path(D4) filesep fname(1:end-4) '.mat'];
    S.to_do=[0 1 1];
    S.modality='EEG';
    S.do_plots=0;
    S.manual_approval=0;
    S.ica_params.approach='defl';
    S.ica_file=[workingdir '/ica_file'];
    fname_out=osl_africa(S);
    D5=spm_eeg_load([fname_out(1:end-4) '.mat']);
    
    
    
    % re-detect bad channels and repair them by interpolation
    
    eeg_list=find(strcmp(chantype(D5),'EEG'));
    bad_list=badchannels(D5);
    variance=var(D5(eeg_list,:,:),[],2);
    bad_channel_var=find(variance>mean(variance(variance>0))+2*std(variance(variance>0)));
    kurt=kurtosis(D5(eeg_list,:,:),0,2);
    bad_channel_kurt=find(kurt>mean(kurt(variance>0))+2*std(kurt(variance>0)));
    list_bch=unique([bad_channel_var ; bad_channel_kurt]);
    vect_bch=zeros(1,size(D5,1));
    vect_bch(bad_list)=1;
    vect_bch(eeg_list(list_bch))=1;
    D5=badchannels(D5,[],vect_bch);
    D5.save;
    
    S=[];
    S.D=D5;
    S.modality='EEG';
    D6=osl_channelrepair(S);
    
    
    % referencing using average reference computing on the good channels only
    
    labelorg = D6.chanlabels;
    good_channel=ones(1,D6.nchannels);
    good_channel(badchannels(D6))=0;
    eeg_list=find(strcmp(chantype(D6),'EEG') & good_channel);
    
    tra = eye(D6.nchannels);
    tra (eeg_list,:)= detrend(tra(eeg_list,:), 'constant');
    labelnew =labelorg;
    
    S = [];
    S.D = D6;
    S.montage.labelorg=labelorg;
    S.montage.labelnew=labelnew;
    S.montage.tra=tra;
    S.keepothers = 'no';
    D7 = spm_eeg_montage(S);
    
    oslview_eeg(D7);
    
    
    % detecting bad events and channels
    start_time=-1000;  % in ms
    stop_time=1000; % in ms
    trig_info=events(D7);
    Fs_new=fsample(D7);
    clear evtimes;
    clear evlabels;
    for k=1:length(trig_info)
        evtimes(k,:)=trig_info(k).time;
        evlabels(k,:)=str2num(trig_info(k).value);
    end
    trl = round(Fs_new*[(evtimes+start_time/1000),(evtimes+stop_time/1000), repmat(start_time/1000,size(evtimes,1),1)]);
    
    S = [];
    S.D = D7;
    S.bc = 0;
    S.epochinfo.trl = trl;
    S.epochinfo.conditionlabels = evlabels;
    S.epochinfo.padding = 0;
    S.modality='EEG';
    S.artifact_amplitude=5; % default settings
    S.false_rate=0.2; % maximum rate of removed triggers
    
    [bad_events,bad_channels]=osl_detect_badevent_old(S);
    
    S2=[];
    S2.D=D7;
    fname=fnamedat(D7);
    S2.outfile=[path(D7) filesep 'S' fname(1:end-4) '.mat'];
    D8=spm_eeg_copy(S2);
    
    ev=events(D8);
    ev(bad_events)=[];
    D8=events(D8,[],[]);  % delete all events
    D8=events(D8,[],ev);
    
    bch=zeros(1,size(D8,1));
    bch(badchannels(D8))=1;
    bch(bad_channels)=1;
    D8 = badchannels(D8, [], bch);
    D8.save;
    
    S=[];
    S.D=D8;
    S.modality='EEG';
    D9=osl_channelrepair(S);
    
    
    % doing epoching
    
    start_time=-100; % in ms
    stop_time=500;  % in ms
    trig_info=events(D9);
    Fs_new=fsample(D9);
    clear evtimes;
    clear evlabels;
    for k=1:length(trig_info)
        evtimes(k,:)=trig_info(k).time;
        evlabels(k,:)=str2num(trig_info(k).value);
    end
    trl = round(Fs_new*[(evtimes+start_time/1000),(evtimes+stop_time/1000), repmat(start_time/1000,size(evtimes,1),1)]);
    
    
    S = [];
    S.D = D9;
    S.bc = 0;
    S.epochinfo.trl = trl;
    S.epochinfo.conditionlabels = evlabels;
    S.epochinfo.padding = 0;
    D_epoched = spm_eeg_epochs(S);
    D_epoched.epochinfo=S.epochinfo;
    D_epoched.save;
    
    
end;




