function [source_recon_sess source_recon_sess_results sess_report] = osl_prepare_oat_source_recon(source_recon_sess, sess_report)

% [source_recon_sess source_recon_results report] = osl_prepare_oat_source_recon(source_recon_sess, report)
% [source_recon_sess source_recon_results] = osl_prepare_oat_source_recon(source_recon_sess )
%
% Called by oat source recon stage to prepare the data before calling the 
% SPM beamforming toolbox
%
% Carries out:
% 1) Copies SPM MEEG object to oat.source_recon.dirname
% 2) Bandpass filtering
% 3) Epoching (optional)
% 4) Establishes valid time windows, trials and channels
% 5) Normalises different modalities
%
% MWW 2014

if nargin>1
    do_report=1;
else
    sess_report=[];
end

source_recon_sess_results=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check epoching in the SPM MEEG objects passed in
only_epoched_data_provided=0;
if ~isempty(source_recon_sess.D_epoched) && isempty(source_recon_sess.D_continuous)
    only_epoched_data_provided=1;
    source_recon_sess.D=source_recon_sess.D_epoched;
elseif ~isempty(source_recon_sess.D_continuous),
    source_recon_sess.D=source_recon_sess.D_continuous;
else
    error('Need to specify source_recon_sess.D_continuous or source_recon_sess.D_epoched.');
end

% do epoching if needed, but only if any epoch info is provided
do_epoching=0;
if ~only_epoched_data_provided
    if ~isempty(source_recon_sess.D_epoched) || ~isempty(source_recon_sess.epochinfo),
        do_epoching=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

disp(['Preparing source recon stage for ' source_recon_sess.D]);
disp(['Will be designated ' source_recon_sess.session_name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% move copy of SPM MEEG object to work on into the working directory
if ~exist(source_recon_sess.dirname,'dir')
    mkdir(source_recon_sess.dirname);
end

% move the spm file(s), giving it a unique name in the parent directory
S2=[];
S2.D=source_recon_sess.D;
tempstring = tempname;
tempstring = tempstring(end-12:end);
S2.outfile=[source_recon_sess.session_name tempstring '_spm_meeg'];
[p spmname e] = fileparts(S2.outfile);                               
S2.outfile=[spmname '.mat'];   

S2.updatehistory=0;
D = spm_eeg_copy(S2);

movefile(fullfile(D.path,D.fname),source_recon_sess.dirname);
movefile(D.fnamedat,source_recon_sess.dirname);
spm_filename=fullfile(source_recon_sess.dirname,D.fname);
D=spm_eeg_load(spm_filename);

% copy the file in the analysis directory in order to give it the standard
% name
S2=[];
S2.D = D;
S2.outfile = [source_recon_sess.session_name '_spm_meeg'];
S2.updatehistory=0;
D = spm_eeg_copy(S2);

% delete the temporary files
matfn = fullfile(source_recon_sess.dirname,[source_recon_sess.session_name tempstring '_spm_meeg.mat']);
datfn = fullfile(source_recon_sess.dirname,[source_recon_sess.session_name tempstring '_spm_meeg.dat']);
delete(matfn);
delete(datfn);

spm_filename = fullfile(source_recon_sess.dirname,[source_recon_sess.session_name '_spm_meeg']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do bandpass filtering
isValidFreqRange = (length(source_recon_sess.freq_range) > 1) && ...
                   (source_recon_sess.freq_range(2) > source_recon_sess.freq_range(1));
if isValidFreqRange
    disp('Temporal filtering...');
    S2=[];
    S2.D=spm_filename;
    S2.freq  = source_recon_sess.freq_range;
    S2.band = 'bandpass';

    % use a zero-lag butterworth iir filter
    S2.type  = 'but';
    S2.order = 5; % 5 is the default in spm_eeg_filter
    S2.dir   = 'twopass';

    Dnew = spm_eeg_filter(S2);   
    
    D.delete;
    D=Dnew;
elseif length(source_recon_sess.freq_range) > 1
    error(['Frequency range incorrectly specified. ', ...
           'Frequencies in source_recon_sess.freq_range: %d must be greater than %d. \n'], ...
          source_recon_sess.freq_range(2), source_recon_sess.freq_range(1));
else
    disp('No temporal filtering applied \n\n');
end%if do bandpass filtering

if source_recon_sess.bandstop_filter_mains
    disp('Applying mains notch filter ....');
    % do notch filtering to remove mains noise
    notchBands = [48 52; 98 102; 148 152; 198 202; 248 252; 298 302; 348 352; 398 402; 448 452; 498 502];
    lowIn      = notchBands(:,1) > source_recon_sess.freq_range(1);
    highIn     = notchBands(:,2) < source_recon_sess.freq_range(2);
    doNotch    = find(lowIn & highIn);
    
    for iNotch = 1:numel(doNotch)
        S3              = [];
        S3.D            = Dnew;
        S3.filter.PHz   = notchBands(doNotch(iNotch),:);
        S3.filter.dir   = 'twopass';
        S3.filter.band  = 'stop';
        disp(['Band-stop filtering from ' num2str(S3.freq(1)) 'Hz to ' num2str(S3.freq(2)) 'Hz.']);
        Dnew = spm_eeg_filter(S3);
    end % for iNotch = 1:numel(doNotch)
    D.delete;
    D=Dnew;

end %if source_recon_sess.bandpass_filter_mains

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do epoching

% epochinfo will be in source_recon_sess.epochinfo, if its empty then look in D_epoched for D_epoched.epochinfo
if(do_epoching)
    disp('Epoching...');
        
    if(~isempty(source_recon_sess.epochinfo))
        S2=source_recon_sess.epochinfo; 
    else
        try
            D_epoched=spm_eeg_load(source_recon_sess.D_epoched);
            S2=D_epoched.epochinfo;
        catch
            error('No epoch info available');
        end
    end
    S2.D = D;
    S2.epochinfo.padding = 0;
    S2.save=0;
    S2.reviewtrials=0;
    S2.bc=0; 
    disp('Doing no within-trial baseline correction at the point of epoching');

    Dnew = spm_eeg_epochs(S2);

    Dnew.save;
    D.delete;
    D=Dnew;
    
    %% get bad channels and trials from passed in source_recon_sess.D_epoched    
    if(~isempty(source_recon_sess.D_epoched))
        D_epoched_passed_in=spm_eeg_load(source_recon_sess.D_epoched);
        rejtmp=badtrials(D_epoched_passed_in);
        rej=zeros(1,D_epoched_passed_in.ntrials);
        rej(rejtmp)=1;
        rej(D_epoched_passed_in.badtrials)=1;
        D = badtrials(D, 1:length(rej), rej);  

        if(length(D_epoched_passed_in.badchannels)>0)
            D = badchannels(D, D_epoched_passed_in.badchannels, ones(length(D_epoched_passed_in.badchannels),1));
        end
        D.save;
    end

end%if do_epoching

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Establish setup channels
if strcmp(source_recon_sess.modalities{1},'EEG')
    modality_meeg='EEG';
else
    modality_meeg='MEGANY';
end
% 
% chanindmeg = strmatch(modality_meeg, D.chantype);
% 
% chanind = setdiff(chanindmeg, D.badchannels);
% if isempty(chanind)
%     error(['No good ' modality_meeg ' channels were found.']);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Establish time windows of interest

if isempty(source_recon_sess.time_range)
    source_recon_sess.time_range = [D.time(1) D.time(end)];
end
time_range = source_recon_sess.time_range;

if(D.ntrials==1)
    goodsamples = good_samples(D,D.indchantype(modality_meeg));
else
    goodsamples = true(1,D.nsamples); 
end

samples_of_interest=zeros(1,D.nsamples);
for i=1:size(time_range,1),
    samples_of_interest(D.indsample(time_range(i, 1)):D.indsample(time_range(i, 2)))=1;
end

samples2use = samples_of_interest & goodsamples;
woi=[D.time(find(diff([0 samples2use])==1))' D.time(find(diff([samples2use 0])==-1))'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Establish trials

if strcmp(source_recon_sess.conditions{1},'all')
    trials = D.indtrial(D.condlist,'good');
else
    
    trials = D.indtrial(source_recon_sess.conditions,'good');
    if isempty(trials)
        error('No trials matched the selection, check the specified condition labels');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalise modalities using, e.g., mean or smallest eigenvalues
%% calculated using good channels and good trials, and over all woi

if strcmp( source_recon_sess.normalise_method,'none')
    S=[];
    S.D=D_epoched;
    S.normalise_method='mean_eig';
    S.datatype='neuromag';
    S.do_plots=true;
    [D_epoched pcadim normalisation]= osl_normalise_sensor_data(S);
else
    chanind=indchantype(D,modality_meeg,'good');
    pcadim=length(chanind);
    normalisation=1;        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Legacy code for HMM

NK=1;

hmm_class=zeros(1,size(D,2),size(D,3));
hmm_class_probs=zeros(1,size(D,2),size(D,3));

hmm_class(1,find(samples2use),trials)=1;
hmm_class_probs(1,find(samples2use),trials)=1;

% add class channels (legacy code)
Sc=[];
Sc.D=D;
Sc.newchandata=[hmm_class; hmm_class_probs];
Sc.newchanlabels{1}='Class';
Sc.newchantype{1}='CLASS';
oldname=D.fullfile;
Sc.newname=prefix(oldname,'concat');
for kk=1:NK
    Sc.newchanlabels{kk+1}=['ClassPr' num2str(kk)];
    Sc.newchantype{kk+1}='CLASSPR';
end

Dnew = osl_concat_spm_eeg_chans( Sc );

D.delete;

% restore file name
Sc=[];
Sc.D=Dnew;
Sc.outfile=oldname;
D=spm_eeg_copy(Sc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

source_recon_sess_results.BF.data.D=spm_eeg_load(D);    
source_recon_sess_results.woi=woi;
source_recon_sess_results.samples2use=samples2use;
source_recon_sess_results.pca_order = pcadim;
source_recon_sess_results.normalisation = normalisation;
source_recon_sess_results.trials=trials;
%source_recon_results.chanind=chanind;

end

