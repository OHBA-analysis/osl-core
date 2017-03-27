% HMM-MAR DEMO SCRIPT FOR INFERRING A GROUP SPECTRALLY DEFINED 
% HMM FROM SOURCE SPACE MEG DATA
% This follows the paper: Vidaurre et al, NeuroImage (2016)
% Diego Vidaurre, Feb 2017
   
% Software directories
HMMMAR_dir = fullfile(osldir,'HMM-MAR');
ohbaexternal_dir = fullfile(osldir,'ohba-external');

% Directory of the data:
data_dir = '/home/diegov/MATLAB/OSL_course/data/HMM-MAR/';
% Name for this HMM-MAR analysis:
hmmmar_name = '/home/diegov/MATLAB/OSL_course/Tutorials/HMMMAR_demo/hmmmar_demo.mat';

% add the HMMMAR paths
% HMM-MAR is standalone, so it does not hold any dependence with OSL or
% other software packages
addpath(genpath(HMMMAR_dir))

% check that PCA is Matlab's own, otherwise remove the corresponding path
which pca

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
    file = [data_dir 'sub' num2str(j) '.mat'];
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
options.covtype = 'full'; % model connectivity
options.zeromean = 0; % model the mean, i.e. model the "amount of power" 
options.cyc = 100;  
options.initcyc = 10;  
options.initrep = 3;  
options.verbose = 1; % show progress?
if do_analysis
    [hmm_env,Gamma_env] = hmmmar(X,T,options);
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
end

%% 1) HMM-MAR on raw signals 
% Such as in Vidaurre et al, NeuroImage (2016)
%% Loading and preparing the data

% We are going to concatenate data into a single matrix
% (if data is too big, the file names can also be provided to the hmmmar function)
X = []; % data, time by regions
T = []; % length of trials, a vector containing how long is each session/trial/subject

subjects = [1 2 4 5 6 7 8 9]; % index of the subjects
N = length(subjects);

for j = subjects % iterate through subjects
    file = [data_dir 'sub' num2str(j) '.mat'];
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
end

%% Because the number of cycles was set to such a low number, 
% we reload a previously computed run, which would have taken a bit longer

load('hmmmar_demo.mat')


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
    file = [data_dir 'sub' num2str(j) '.mat'];
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
    file = [data_dir 'sub' num2str(j) '.mat'];
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
subplot(1,2,1)
plot((-halfwindow:halfwindow)/Hz,evokedGamma_env,'LineWidth',2)
xlabel('Time (s)','FontSize',15)
ylabel('State probability','FontSize',15)
title('HMM on envelope','FontSize',17)
ylim([0.05 0.7])
subplot(1,2,2)
plot((-halfwindow:halfwindow)/Hz,evokedGamma_raw,'LineWidth',2)
xlabel('Time (s)','FontSize',15)
ylabel('State probability','FontSize',15)
title('HMM-MAR on raw signals','FontSize',17)
ylim([0.05 0.7])

%% Show the spectral info (power and coherence) for the MAR

colors = {'b',[0.2 0.7 0.2],'r'};

figure(2);clf(2) % parametric spectra, from the HMM-MAR run
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


figure(4);clf(4) % parametric spectra, from the HMM-MAR run
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

[psd_tf,coh_tf] = hmmtimefreq(spectra_raw,evokedGamma_raw,1);
f = spectra_raw.state(1).f ;
indf = f>4 & f<=30;
f = f(indf);

figure(5);clf(5)
subplot(3,1,1)
l = max(max(max(abs(psd_tf(indf,:,1)))),max(max(abs(psd_tf(indf,:,1)))));
imagesc((-halfwindow:halfwindow)/Hz,f,psd_tf(:,indf,1)',[-l l]);colorbar
hold on; plot([0 0],[4 30],'k'); hold off
xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
title('Power channel 1','FontSize',16)
subplot(3,1,2)
imagesc((-halfwindow:halfwindow)/Hz,f,psd_tf(:,indf,2)',[-l l]);colorbar
hold on; plot([0 0],[4 30],'k'); hold off
xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
title('Power channel 2','FontSize',16)
subplot(3,1,3)
l = max(max(abs(coh_tf(indf,:,1,2))));
imagesc((-halfwindow:halfwindow)/Hz,f,coh_tf(:,indf,1,2)',[-l l]);colorbar
hold on; plot([0 0],[4 30],'k'); hold off
xlabel('Time (s)','FontSize',14); ylabel('Frequency (Hz)','FontSize',14);
title('Coherence','FontSize',16)


%% Probability of transit between states

P = hmm_raw.P; 
for k=1:3
    P(k,k) = 0;
    P(k,:) = P(k,:) / sum(P(k,:));
end
P

figure(6)
imagesc(P);colorbar
title('Transition probability matrix','FontSize',16)
set(gca,'Xtick',1:3,'ytick',1:3,'FontSize',14)
xlabel('To','FontSize',16)
ylabel('From','FontSize',16)