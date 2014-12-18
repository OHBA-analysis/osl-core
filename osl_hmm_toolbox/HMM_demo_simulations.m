% Infer covariance structure and state time courses for simulated data with
% three distinct states using a hidden Markov model


%% HMM Simulation 1: variance changes

% state 1 variance increased in 1: 
Sigma = [2 0 0 0 
         0 1 0 0
         0 0 1 0
         0 0 0 1];
k1 = randn(1000,4)*chol(Sigma);

% state 2 variance increased in 2: 
Sigma = [1 0 0 0
         0 2 0 0
         0 0 1 0
         0 0 0 1];
k2 = randn(1000,4)*chol(Sigma);

% state 3 variance increased in 3: 
Sigma = [1 0 0 0
         0 1 0 0
         0 0 2 0
         0 0 0 1];
k3 = randn(1000,4)*chol(Sigma);

% concatenate temporally:
simdata = [k1;k2;k3];



% infer HMM with 3 states:
options.K      = 3; % number of states
options.order  = 0; % set 0 to run Multivariate normal observation model
options.Ninits = 1; % Number of initialisations
options.Hz     = 1; % sampling rate
[hmm,stats] = osl_hmm_infer(simdata,options);

% Plot simulated data
figure, subplot(4,3,1:3), plot(simdata), axis tight; set(gca,'Title',text('String','Simulated data'))

% Plot simulated covariance matrices
subplot(4,3,4), imagesc(cov(k1)), colormap('gray')
subplot(4,3,5), imagesc(cov(k2)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(4,3,6), imagesc(cov(k3)), colormap('gray')

% Plot inferred covariance matrices
subplot(4,3,7), imagesc(hmm.state(1).Cov), colormap('gray')
subplot(4,3,8), imagesc(hmm.state(2).Cov), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(4,3,9), imagesc(hmm.state(3).Cov), colormap('gray')

% Plot state path
subplot(4,3,10:12); hold on
plot(hmm.statepath), set(gca,'Title',text('String','Inferred state path'))
set(gca,'ylim',[0.5 3.5]); ylabel('state #')

set(gcf,'pos',[100 100 600 800],'color','w')




%% HMM Simulation 2: covariance changes


% state 1 covariance increased between 1 & 2: 
Sigma = [1    0.5  0   0
         0.5  1    0   0
         0    0    1   0
         0    0    0   1];
k1 = randn(1000,4)*chol(Sigma);

% state 2 covariance increased between 1 & 3: 
Sigma = [1    0   0.5  0
         0    1   0    0
         0.5  0   1    0
         0    0   0    1];
k2 = randn(1000,4)*chol(Sigma);

% state 3 covariance increased between 2 & 4: 
Sigma = [1   0     0    0
         0   1     0    0.5
         0   0     1    0
         0   0.5   0    1];
k3 = randn(1000,4)*chol(Sigma);

% concatenate temporally:
simdata = [k1;k2;k3];


% infer HMM with 3 states:
options.K      = 3; % number of states
options.order  = 0; % set 0 to run Multivariate normal observation model
options.Ninits = 1; % Number of initialisations
options.Fs     = 1; % sampling rate
[hmm,stats] = osl_hmm_infer(simdata,options);

% Plot simulated data
figure, subplot(4,3,1:3), plot(simdata), axis tight; set(gca,'Title',text('String','Simulated data'))

% Plot simulated covariance matrices
subplot(4,3,4), imagesc(cov(k1)), colormap('gray')
subplot(4,3,5), imagesc(cov(k2)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(4,3,6), imagesc(cov(k3)), colormap('gray')

% Plot inferred covariance matrices
subplot(4,3,7), imagesc(hmm.state(1).Cov), colormap('gray')
subplot(4,3,8), imagesc(hmm.state(2).Cov), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(4,3,9), imagesc(hmm.state(3).Cov), colormap('gray')

% Plot state path
subplot(4,3,10:12); hold on
plot(hmm.statepath), set(gca,'Title',text('String','Inferred state path'))
set(gca,'ylim',[0.5 3.5]); ylabel('state #')

set(gcf,'pos',[100 100 600 800],'color','w')