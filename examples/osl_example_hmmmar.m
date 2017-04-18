%% HMM-MAR 

% This example shows how to use the HMM-MAR to infer transient states that:
% (i) are pectrally defined, i.e. the characteristics of interest are defined as a function of frequency.
% (ii) are based on the raw time series, i.e. we do not need to bandpass filter or compute power envelopes.
% (iii) are not only sensitive to power differences but also to phase coupling.

% The script infers a group (spectrally-defined) HMM from source space MEG data, following the paper
% Vidaurre et al, NeuroImage (2016)
   
% Directory of the data:
data_dir = fullfile(osldir,'example_data','hmmmar_example');
% Name for this HMM-MAR analysis:
hmmmar_name = fullfile(osldir,'example_data','hmmmar_example','hmmmar_demo.mat');

% Set |do_analysis=1| to re-run the analysis, otherwise use precomputed result
do_analysis = 0; 

%% 1) HMM-MAR on hilbert envelopes 
% Such as in Baker et al, Elife (2014)
%% Loading and preparing the data

% We are going to concatenate data into a single matrix
% (if data is too big, the file names can also be provided to the hmmmar function)
X = []; % data, time by regions
T = []; % length of trials, a vector containing how long is each session/trial/subject

subjects = [1 2 4 5 6 7 8 9]; % index of the subjects
N = length(subjects);

for j = subjects % iterate through subjects
    file = fullfile(data_dir,['sub' num2str(j) '.mat']);
    load(file); % load the data
    sourcedata = sourcedata' ; % we need it (time by channels)
    T = [T size(sourcedata,1)]; % time length of this subject
    for i=1:size(sourcedata,2)
        sourcedata(:,i) = abs(hilbert(sourcedata(:,i))); % Hilbert envelope
    end
    sourcedata = zscore(sourcedata); % standardize subject such that it has mean 0 and stdev 1
    X = [X; sourcedata]; % concat subject
end

%% Prepare the HMM options and run the HMMMAR
    
options = struct();
options.K = 3; % number of states
options.order =  0; % order 0 means a Gaussian distribution (adequate for power time series)
options.covtype = 'full'; 
% model connectivity, covtype='full' means that we model the covariance between regions
options.zeromean = 0; 
% zeromean=0 means that we model the mean, i.e. model the "amount of power"
% zeromean=1 means that we do *not* model the mean
options.cyc = 100;  
options.initcyc = 10;  
options.initrep = 3;  
options.verbose = 1; % show progress?
if do_analysis
    [hmm_env,Gamma_env] = hmmmar(X,T,options);
else
    load(hmmmar_name,'hmm_env','Gamma_env')
end

%% Compute the spectral information of the states (power, coherence, etc) 
% using a weighted version of the multitaper
% Note that we can't get an estimation of the spectra directly from the
% parameters in this case, as the HMM has been run in wideband power time series

options = struct();
options.fpass = [1 40]; % frequency range we want to look at
options.tapers = [4 7]; % internal multitaper parameter
options.Fs = 200; % sampling frequency in hertzs
options.win = 10 * options.Fs; % window, related to the level of detail in frequency of the estimation

if do_analysis
    spectra_env = hmmspectramt(X,T,Gamma_env,options);
else
    load(hmmmar_name,'spectra_env')
end

%% 2) HMM-MAR on raw signals 
% Such as in Vidaurre et al, NeuroImage (2016)
%% Loading and preparing the data

% We are going to concatenate data into a single matrix
% (if data is too big, the file names can also be provided to the hmmmar function)
X = []; % data, time by regions
T = []; % length of trials, a vector containing how long is each session/trial/subject

subjects = [1 2 4 5 6 7 8 9]; % index of the subjects
N = length(subjects);

for j = subjects % iterate through subjects
    file = fullfile(data_dir,['sub' num2str(j) '.mat']);
    load(file); % load the data
    sourcedata = sourcedata' ; % we need it (time by channels)
    T = [T size(sourcedata,1)]; % time length of this subject
    sourcedata = zscore(sourcedata); % standardize subject such that it has mean 0 and stdev 1
    X = [X; sourcedata]; % concat subject
end

%% Signs might be arbitrarily flipped, so we need to disambiguate this
% this is because the intrinsic ambiguity in source-reconstructed signals

options = struct();
options.maxlag = 8;
options.noruns = 100;

[flips,scorepath] = findflip(X,T,options);
X = flipdata(X,T,flips);

%% And run the HMM-MAR on this

options = struct();
options.K = 3; % number of states
options.order =  4; % MAR order 
options.covtype = 'full';
options.cyc = 100;  
options.initcyc = 10;  
options.initrep = 3;
options.verbose = 1; % show progress?
if do_analysis
    [hmm_raw,Gamma_raw] = hmmmar(X,T,options);
else
    load(hmmmar_name,'hmm_raw','Gamma_raw')
end

%% Compute the spectral information of the states (power, coherence, etc) 
% Because a MAR model implicitly contains all the spectral information, 
% we can use the parameters to derive power, coherence, phase relations, etc.

options = struct();
options.fpass = [1 40]; % frequency range we want to look at
options.Nf = 100; % number of frequency bins
options.Fs = 200; % sampling frequency in hertzs
options.order = 11;

if do_analysis
    spectra_raw = hmmspectramar(X,T,[],Gamma_raw,options);
else
    load(hmmmar_name,'spectra_raw')
end

%% Now we do some plotting of the HMM on power envelopes and HMM-MAR results
% computed above

% Because the number of cycles was set to such a low number, 
% we reload a pre-computed run, which would have taken a bit longer

load(hmmmar_name)

%% Now let's see the state evoked probability, locked to the stimulus
% the stimulus is saved in the variable onset, which contains a (Tx1)
% logical vector with 1 when the fingertapping is effected. 

subjects = [1 2 4 5 6 7 8 9]; % index of the subjects
Hz = 200; % sampling frequency
L = 8; % length of the window around the stimulus (seconds)
window = Hz*L+1;

% "evoked state response" around the stimulus, for the envelope run
evokedGamma_env = zeros(window,hmm_env.K,length(subjects));
t0 = 0;
for j = subjects % iterate through subjects
    file = fullfile(data_dir,['sub' num2str(j) '.mat']);
    load(file,'onset'); % load the data
    T = length(onset);
    index = t0 + (1:T);
    Gamma_subj = Gamma_env(index,:); % state time course of this subject
    evokedGamma_env(:,:,j) = evokedStateProbability(onset,T,Gamma_subj,window); 
    t0 = t0 + length(index);
end
evokedGamma_env = mean(evokedGamma_env,3); % average across subjects

% "evoked state response" around the stimulus, for the envelope run
evokedGamma_raw = zeros(window,hmm_env.K,length(subjects));
t0 = 0;
for j = subjects % iterate through subjects
    file = fullfile(data_dir,['sub' num2str(j) '.mat']);
    load(file,'onset'); % load the data
    T = length(onset);
    index = t0 + (1:T-hmm_raw.train.order);
    Gamma_subj = Gamma_raw(index,:); % state time course of this subject
    evokedGamma_raw(:,:,j) = evokedStateProbability(onset,T,Gamma_subj,window);
    t0 = t0 + length(index);
end
evokedGamma_raw = mean(evokedGamma_raw,3); % average across subjects
 
% And plot both 
figure(1);
halfwindow = (window-1)/2;
% Evoked state response for the HMM on power envelopes
subplot(1,2,1)
plot((-halfwindow:halfwindow)/Hz,evokedGamma_env,'LineWidth',2)
xlabel('Time (s)','FontSize',15)
ylabel('State probability','FontSize',15)
title('HMM on envelope','FontSize',17)
ylim([0.05 0.7])
% Evoked state response for the HMM-MAR
subplot(1,2,2)
plot((-halfwindow:halfwindow)/Hz,evokedGamma_raw,'LineWidth',2)
xlabel('Time (s)','FontSize',15)
ylabel('State probability','FontSize',15)
title('HMM-MAR on raw signals','FontSize',17)
ylim([0.05 0.7])

%% Show the spectral info (power and coherence) for the 
% HMM on power envelopes and HMM-MAR

colors = {'b',[0.2 0.7 0.2],'r'};

% Non-parametric spectra from the HMM on power envelopes
figure(2);clf(2)  
for k=1:hmm_env.K
    subplot(2,2,1)
    hold on
    plot(spectra_env.state(k).f,spectra_env.state(k).psd(:,1,1),'Color',colors{k},'LineWidth',2);xlim([4 30])
    set(gca,'FontSize',14)
    hold off
    subplot(2,2,4)
    hold on
    plot(spectra_env.state(k).f,spectra_env.state(k).psd(:,2,2),'Color',colors{k}','LineWidth',2);xlim([4 30])
    set(gca,'FontSize',14)
    hold off
    subplot(2,2,3)
    hold on
    plot(spectra_env.state(k).f,spectra_env.state(k).coh(:,1,2),'Color',colors{k},'LineWidth',2);xlim([4 30])
    set(gca,'FontSize',14)
    hold off
end
subplot(2,2,1)
title('Power Chan1, HMM-Gaussian','FontSize',16)
subplot(2,2,4)
title('Power Chan2, HMM-Gaussian','FontSize',16)
subplot(2,2,3)
title('Coherence, HMM-Gaussian','FontSize',16)

% Parametric spectra (MAR-based) from the HMM-MAR
figure(3);clf(3)  
for k=1:hmm_env.K
    subplot(2,2,1)
    hold on
    plot(spectra_raw.state(k).f,spectra_raw.state(k).psd(:,1,1),'Color',colors{k},'LineWidth',2);xlim([4 30])
    set(gca,'FontSize',14)
    hold off
    subplot(2,2,4)
    hold on
    plot(spectra_raw.state(k).f,spectra_raw.state(k).psd(:,2,2),'Color',colors{k}','LineWidth',2);xlim([4 30])
    set(gca,'FontSize',14)
    hold off
    subplot(2,2,3)
    hold on
    plot(spectra_raw.state(k).f,spectra_raw.state(k).coh(:,1,2),'Color',colors{k},'LineWidth',2);xlim([4 30])
    set(gca,'FontSize',14)
    hold off
end
set(gca,'FontSize',14)
subplot(2,2,1)
title('Power Chan1, HMM-MAR','FontSize',16)
subplot(2,2,4)
title('Power Chan2, HMM-MAR','FontSize',16)
subplot(2,2,3)
title('Coherence, HMM-MAR','FontSize',16)


%% HMM-based Time-frequency analysis
% Given the higher sensitivity of the HMM-MAR spectra, apparent in the
% previous block, we focus here just on the HMM-MAR results 

[psd_tf,coh_tf] = hmmtimefreq(spectra_raw,evokedGamma_raw,1);
f = spectra_raw.state(1).f ;
indf = f>4 & f<=30;
f = f(indf);

figure(4);clf(4)
subplot(3,1,1)
l = max(max(max(abs(psd_tf(indf,:,1)))),max(max(abs(psd_tf(indf,:,1)))));
imagesc((-halfwindow:halfwindow)/Hz,f,psd_tf(:,indf,1)',[-l l]);colorbar
hold on; plot([0 0],[4 30],'k'); hold off
xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
title('Power channel 1 : HMM-MAR','FontSize',16)
subplot(3,1,2)
imagesc((-halfwindow:halfwindow)/Hz,f,psd_tf(:,indf,2)',[-l l]);colorbar
hold on; plot([0 0],[4 30],'k'); hold off
xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
title('Power channel 2 : HMM-MAR','FontSize',16)
subplot(3,1,3)
l = max(max(abs(coh_tf(indf,:,1,2))));
imagesc((-halfwindow:halfwindow)/Hz,f,coh_tf(:,indf,1,2)',[-l l]);colorbar
hold on; plot([0 0],[4 30],'k'); hold off
xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
title('Coherence : HMM-MAR','FontSize',16)


%% Probability of transit between states for the 
% HMM on power envelopes and HMM-MAR

P_env = hmm_env.P; 
% We normalize, such that we focus on the between state transitions
for k=1:3
    P_env(k,k) = 0;
    P_env(k,:) = P_env(k,:) / sum(P_env(k,:));
end
disp('HMM-Gaussian Probability of transition from state i to state j')
P_env


P_raw = hmm_raw.P; 
% We normalize, such that we focus on the between state transitions
for k=1:3
    P_raw(k,k) = 0;
    P_raw(k,:) = P_raw(k,:) / sum(P_raw(k,:));
end
disp('HMM-MAR Probability of transition from state i to state j')
P_raw



figure(5)

% plot trans prob matrix for the HMM-Gaussian
subplot(1,2,1)
imagesc(P_env);colorbar
title('Transition probability matrix for the HMM-Gaussian','FontSize',16)
set(gca,'Xtick',1:3,'ytick',1:3,'FontSize',14)
xlabel('To','FontSize',16)
ylabel('From','FontSize',16)

% plot trans prob matrix for the HMM-MAR
subplot(1,2,2)
imagesc(P_raw);colorbar
title('Transition probability matrix for the HMM-MAR','FontSize',16)
set(gca,'Xtick',1:3,'ytick',1:3,'FontSize',14)
xlabel('To','FontSize',16)
ylabel('From','FontSize',16)


