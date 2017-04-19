%%
% This example shows how to use the HMM-MAR to infer transient states that:
%
% (i) are spectrally defined, i.e. the characteristics of interest are defined as a function of frequency,
% (ii) are based on the raw time series, i.e. we do not need to bandpass filter or compute power envelopes,
% (iii) are not only sensitive to power differences but also to phase coupling.
%
% In particular, this script will estimate a group (spectrally-defined) HMM from MEG task data, 
% using two source-reconstructed regions in the motor cortex (left and right).
% The task is a finger-tapping motor task, where subjects press a button volitionally.
% We will see whether, with no prior information of the task during training, 
% the estimated HMM contains states that are temporally related to the task timing,
% and we will inspect if these states contain meaningful features.
% We will compare the HMM-MAR estimation with the estimation of an HMM on the power time series (Baker et al, 2014).
%
% The script follows the paper Vidaurre et al. (2016)

% Directory of the data:
data_dir = fullfile(osldir,'example_data','hmmmar_example');
% Name for this HMM-MAR analysis:
hmmmar_name = fullfile(osldir,'example_data','hmmmar_example','hmmmar_demo.mat');
% Set do_analysis=1 to re-run the analysis, otherwise use precomputed result
do_analysis = 0; 

%% HMM-Gaussian on hilbert envelopes 
% First we are running the HMM on the Hilbert envelopes (or power time series) using a
% Gaussian distribution as observation model.
% This is the method established in Baker et al. (2014), but applied in task and in two regions only.

%% 
% We first need to load and prepare the data.
% Data, for any variety of the HMM, is usually provided in a matrix (time by regions) with all concatenated subjects.
% However, if data is too big, we can pass a cell (subjects by 1) where each elements contains
% the name of the file containing the data for one subject. 

% Here, we are going to concatenate data into a single matrix. 
% We will also need to let the toolbox know the length of each subject's data (within this big matrix).
X = []; % data, time by regions
T = []; % length of trials, a vector containing how long is each session/trial/subject

subjects = [1 2 4 5 6 7 8 9]; % index of the subjects that we are using
N = length(subjects);

for j = subjects % iterate through subjects
    file = fullfile(data_dir,['sub' num2str(j) '.mat']);
    load(file); % load the data
    sourcedata = sourcedata' ; % data is stored as (channels by time); we need it (time by channels)
    T = [T size(sourcedata,1)]; % time length of this subject
    X = [X; sourcedata]; % concat subject
end

%% 
% Prepare the HMM options, where we will configure it to use  and run the HMMMAR
    
options = struct();
options.K = 3; 
% number of states: in general, the higher this number, the more "detailed" will be the segmentation.
% (By running the model with different number of states, we could get some sense of "state hierarchy").
options.order =  0; 
% order=0 corresponds to a Gaussian distribution (adequate for power time series, i.e. Baker et al's approach).
% By setting order>0, we will be running the HMM-MAR (see below).
options.covtype = 'full'; 
% model connectivity, covtype='full' means that we model the covariance between regions;
% otherwise, we would set covtype='diag', such that we would ignore the convariance and would focus on variance.
options.zeromean = 0; 
% zeromean=0 means that we model the mean, i.e. model the "amount of power"; 
% zeromean=1 means that we do *not* model the mean.
options.standardise = 1; 
% standardize each subject such that it has mean 0 and stdev 1;
% this is typically done in order to avoid between-subject differences being the main driving cause of the states.
options.onpower = 1; 
% onpower=1 indicates that we will be taking the Hilbert envelopes of the signal (ignoring phase)
options.cyc = 100; 
% number of training cycles, although might stop earlier if convergence is attained. 
options.verbose = 1; 
% show progress?

%% 
% And run the HMM-MAR. Note that the same function will be used for both "HMM-flavours".
if do_analysis % run the HMM-MAR
    [hmm_env,Gamma_env] = hmmmar(X,T,options);
else % load a precomputed run
    load(hmmmar_name,'hmm_env','Gamma_env')
end

%% 
% Now that we have the HMM model, we shall compute the spectral information of the states (power, coherence, etc).
% In the case of the HMM-MAR, as we will see below, the MAR parameters themselves contain this spectral information.
% Here, because we are using a Gaussian distribution on wideband power, 
% we don't know what is happening in the frequency
% domain when each state is active, above and beyond gross changes in wideband power. 
% For example, in our current HMM estimation there is no information about phase coupling.
% Hence, we need to estimate this frequency information (including phase coupling), 
% looking back at the data and using the state time courses we inferred. 
% For this, we will be using a weighted version of the multitaper, proposed in (Vidaurre et al. 2016).
% In brief, this basically weights the data with the state time course for each state, 
% such that segments of the data where a given state has a higher probability of being active will 
% contribute more to the final spectral estimation.

% We set the options for the spectral estimation:
options = struct();
options.fpass = [1 40]; % frequency range we want to look at, in this case between 1 and 40 Hertzs.
options.tapers = [4 7]; % internal multitaper parameter
options.Fs = 200; % sampling frequency in Hertzs
options.win = 10 * options.Fs; % window length, related to the level of detail of the estimation;
% that is, if we increase the win parameter, we will obtain an estimation that is more detailed in the frequency scale
% (i.e. contains more frequency bins between 1 and 40 Hertzs) at the expense of some robustness,.

if do_analysis % Estimate the spectra 
    spectra_env = hmmspectramt(X,T,Gamma_env,options);
else % load a precomputed results
    load(hmmmar_name,'spectra_env')
end


%% HMM-MAR on raw time series

% We have estimated a HMM-Gaussian on the power time series, 
% now we look at the raw data (which contains phase) following (Vidaurre et al. 2016).

%% 
% One thing we must take care of when working on raw data in source-space is sign ambiguity.
% As a consequence of the nature of the source-reconstruction process, the sign of the time series 
% are arbitrary. That means that, because different subjects might have different signs in different channels,
% the phase relations between regions can cancel across subjects. 
% The HMM-MAR toolbox includes a set of functions to correct for sign ambiguity. 

options = struct();
options.maxlag = 8;
options.noruns = 100;
options.verbose = 0;

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


