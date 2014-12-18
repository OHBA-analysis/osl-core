function [Gamma,Gammasum,Xi]=hsinference(data,T,hmm,residuals)
%
% inference engine for HMMs.
%
% INPUT
%
% data      Observations - a struct with X (time series) and C (classes)
% T         Number of time points for each time series
% hmm       hmm data structure
% residuals in case we train on residuals, the value of those.
%
% OUTPUT
%
% Gamma     Probability of hidden state given the data
% Gammasum  sum of Gamma over t
% Xi        joint Prob. of child and parent states given the data
%
% Author: Diego Vidaurre, OHBA, University of Oxford


N = length(T);
if hmm.train.order > 0, order = 1:hmm.train.timelag:hmm.train.order; order = order(end); 
else order = 0; end 

K=hmm.K;
Gamma=[];
Gammasum=zeros(N,K);
Xi=[];

for in=1:N
    t0 = sum(T(1:in-1)); s0 = t0 - order*(in-1);
    X = data.X(t0+1:t0+T(in),:);
    C = data.C(t0+1:t0+T(in),:);
    R = [zeros(order,size(residuals,2));  residuals(s0+1:s0+T(in)-order,:)];
    % we jump over the fixed parts of the chain
    t = order+1; 
    xi = []; gamma = []; gammasum = zeros(1,K);
    while t<=T(in)
        if isnan(C(t,1)), no_c = find(~isnan(C(t:T(in),1)));
        else no_c = find(isnan(C(t:T(in),1)));
        end
        if t>order+1
            if isempty(no_c), slice = (t-order-1):T(in); slicer = (t-1):T(in);
            else slice = (t-order-1):(no_c(1)+t-2); slicer = (t-1):(no_c(1)+t-2);
            end;
        else
            if isempty(no_c), slice = (t-order):T(in); slicer = t:T(in);
            else slice = (t-order):(no_c(1)+t-2); slicer = t:(no_c(1)+t-2);
            end;
        end
        if isnan(C(t,1))
            [gammat,xit]=nodecluster(X(slice,:),K,hmm,R(slicer,:));
        else            
            gammat = zeros(length(slicer),K); 
            if t==order+1, gammat(1,:) = C(slicer(1),:); end
            xit = zeros(length(slicer)-1, K^2);
            for i=2:length(slicer)
                gammat(i,:) = C(slicer(i),:);
                xitr = gammat(i-1,:)' * gammat(i,:) ;
                xit(i-1,:) = xitr(:)';
            end
        end
        if t>order+1, gammat = gammat(2:end,:); end
        xi = [xi; xit];
        gamma = [gamma; gammat];
        gammasum = gammasum + sum(gammasum);
        if isempty(no_c), break;
        else t = no_c(1)+t-1; 
        end;
    end
    Gamma=[Gamma;gamma];
    Gammasum(in,:)=gammasum;
    Xi=cat(1,Xi,reshape(xi,T(in)-order-1,K,K));
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gamma,Xi]=nodecluster(X,K,hmm,residuals)
% inference using normal foward backward propagation

T = size(X,1);
P=hmm.P;
Pi=hmm.Pi;
if hmm.train.order > 0, order = 1:hmm.train.timelag:hmm.train.order; order = order(end); 
else order = 0; end 

B = obslike(X,hmm,residuals);

scale=zeros(T,1);
alpha=zeros(T,K);
beta=zeros(T,K);

alpha(1+order,:)=Pi.*B(1+order,:);
scale(1+order)=sum(alpha(1+order,:));
alpha(1+order,:)=alpha(1+order,:)/scale(1+order);
for i=2+order:T
    alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
    scale(i)=sum(alpha(i,:));		% P(X_i | X_1 ... X_{i-1})
    alpha(i,:)=alpha(i,:)/scale(i);
end;

beta(T,:)=ones(1,K)/scale(T);
for i=T-1:-1:1+order
    beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i);
end;
Gamma=(alpha.*beta);
Gamma=Gamma(1+order:T,:);
Gamma=rdiv(Gamma,rsum(Gamma));

Xi=zeros(T-1-order,K*K);
for i=1+order:T-1
    t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
    Xi(i-order,:)=t(:)'/sum(t(:));
end
