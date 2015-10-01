%% Synthetic data

clear 
addpath(genpath('~/Work/Matlab/HMM-MAR'))
% Specify parameters
K = 4; % number of states
T = 10000; % number of data points
N = 3; % number of channels
epsilon = 0.1; %level of noise
StatePermanency = 100; % factor for the diagonal of the transtion matrix
%%
hmmtrue = struct();
hmmtrue.K = K;
hmmtrue.state = struct();
hmmtrue.train.covtype = 'full';
hmmtrue.train.zeromean = 0;
hmmtrue.train.maxorder = 0;
hmmtrue.train.order = 0;
hmmtrue.train.orderoffset = 0;
hmmtrue.train.timelag = 0;
hmmtrue.train.exptimelag = 0;
hmmtrue.train.S =ones(N);
hmmtrue.train.Sind = ones(1,N);
hmmtrue.train.multipleConf = 0;


for k = 1:K
    hmmtrue.state(k).W.Mu_W = rand(1,N);
    hmmtrue.state(k).Omega.Gam_rate = randn(N) + eye(N);
    hmmtrue.state(k).Omega.Gam_rate = epsilon * 1000 * ...
        hmmtrue.state(k).Omega.Gam_rate' * hmmtrue.state(k).Omega.Gam_rate;
    hmmtrue.state(k).Omega.Gam_shape = 1000;
end

hmmtrue.P = ones(K,K) + StatePermanency * eye(K); %rand(K,K);
for j=1:K,
    hmmtrue.P(j,:)=hmmtrue.P(j,:)./sum(hmmtrue.P(j,:));
end;
hmmtrue.Pi = ones(1,K); %rand(1,K);
hmmtrue.Pi=hmmtrue.Pi./sum(hmmtrue.Pi);

[X,T,Gammatrue] = simhmmmar(T,hmmtrue,[]);

%% Train  models

KK = 2:8;
DD = [2 10 20 100 200];
repetitions = 10;

options = struct();
options.covtype = 'full';
options.order = 0;
options.tol = 1e-7;
options.cyc = 100;
options.zeromean = 0;
options.inittype = 'GMM';
options.verbose = 0;

FE = zeros(length(KK),length(DD),repetitions);
Ktrained = zeros(length(KK),length(DD),repetitions);
fefinal = Inf;

for r = 1:repetitions
    for idd=1:length(DD)
        dd = DD(idd);
        for ik=1:length(KK)
            k = KK(ik);
            fprintf('Rep %d; Using K = %d, DD=%d \n',r,k,dd);
            options.K = k;
            options.DirichletDiag = dd;
            [hmm, ~,~, ~, ~, ~, fehist] = hmmmar(X,T,options);
            FE(ik,idd,r) = fehist(end);
            Ktrained(ik,idd,r) = length(hmm.state);
            %if FE(ik,idd)<fefinal
            %    fefinal = FE(ik,idd);
            %    hmmfinal = hmm;
            %    Gammafinal = Gamma;
            %    Kfinal = k;
            %    DDfinal = dd;
            %end
        end
    end
    %fprintf('Selected K = %d; DD=%d \n',Kfinal,DDfinal);
end

% And final run with the selected parameters
meanFE = mean(FE,3);
meanKtrained = mean(Ktrained,3);
[~,I] = min(meanFE(:));
[I1,I2] = ind2sub(size(meanFE),I);
options.K = KK(I1);
options.DirichletDiag = DD(I2);
[hmm, Gamma,~, ~, ~, ~, fehist] = hmmmar(X,T,options);

save('example.mat')


%% Plot results

addpath('../HMM-MAR-scripts/distributionPlot/')
load('example.mat')

figure(1)

subplot(1,2,1); 
imagesc(meanFE); colorbar; 
set(gca,'XTick',1:length(DD))
set(gca,'XTickLabel',DD,'FontSize',16)
set(gca,'YTick',1:length(KK))
set(gca,'YTickLabel',KK,'FontSize',16)
xlabel('Prior on diag(P)','FontSize',18)
ylabel('Initial K','FontSize',18)
title('Free Energy','FontSize',20)

subplot(1,2,2); 
imagesc(meanKtrained); colorbar; 
set(gca,'XTick',1:length(DD))
set(gca,'XTickLabel',DD,'FontSize',16)
set(gca,'YTick',1:length(KK))
set(gca,'YTickLabel',KK,'FontSize',16)
xlabel('Prior on diag(P)','FontSize',18)
ylabel('Initial K','FontSize',18)
title('Final K','FontSize',20)

Colors = {'b','g','r','m','k'};

figure(2); clf(2)
to_show = [1 2 3 4 5];
for d=to_show
    subplot(1,length(to_show),find(d==to_show));
    %boxplot(squeeze(FE(:,d,:))')
    hold on; 
    %distributionPlot(squeeze(FE(:,d,:))','xValues',KK);
    %for r=1:repetitions
    %    plot(KK,FE(:,d,r)','b','LineWidth',1)
    %end
    plot(KK,meanFE(:,d),'b','LineWidth',3)
    errorbar(KK,meanFE(:,d),std(squeeze(FE(:,d,:)),1,2) )
    xlim([KK(1) KK(end)])
    hold off
    title(strcat('dd=',num2str(DD(d))),'FontSize',18)
    xlabel('No. of states','FontSize',18)
    ylabel('FE','FontSize',18)
end
%xlabel('FE','FontSize',18)
%ylabel('K','FontSize',18)
%title('Free energy for the chosen dd','FontSize',20)

figure(3); clf(3)
subplot(2,1,1); 
plot(1:200,Gammatrue(101:300,:)); ylim([-.25 1.25])
xlabel('Time','FontSize',18)
ylabel('True State time courses','FontSize',18)
subplot(2,1,2); 
plot(1:200,Gamma(101:300,:)); ylim([-.25 1.25])
xlabel('Time','FontSize',18)
ylabel('Estimated state time courses','FontSize',18)
