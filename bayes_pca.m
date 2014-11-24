function [m_w, alpha, C, taus, tau, m_x, Sigma_w, Sigma_x, Q, N, D,alphamin,alphamax]=bayes_pca(x,Q,niters,dofs_lost)

% [m_w, alpha, C, taus]=bayes_pca(x,Q,niters,dofs_lost)
% [m_w, alpha, C, taus]=bayes_pca(x)
%
% input x needs to d*Nrepeats
% Q is max number of components
% niters is number of iterations
% dofs_lost is DOF used up to estimate means in the data
%
% output C will be a d*d covariance matrix
% m_w is the d*Q principle components
% alpha is the Q ARD precisions
% taus is a vector of the noise precision at each iteration
%
% E.g:
% T=1000;N=200;alpha=2/N*randn(N,N);alpha=alpha*alpha';alpha=alpha-diag(diag(alpha));x=mvnrnd(zeros(N,1),inv(eye(N)-alpha),T)';cov(x');[m_w, alpha, C]=bayes_pca2(x,100);

if(nargin<3)
    niters=30;
    dofs_lost=1; % for estimating the mean
end;

if(0)
    T=1000;N=200;RNK=50;
    y=randn(RNK,T);S=randn(N,RNK);
    x=(S*y)';
    [m_w, alpha, C, taus, tau, m_x, Sigma_w, Sigma_x, Q, N, D,alphamin,alphamax]=bayes_pca(x,100);
    figure;plot(alpha);
    min(find(alpha>sqrt(alphamax*alphamin))-1)
end;

% demean and normalise:
t=demean(x,2);
nor=std(squash(t));
t=t/nor;

N=size(t,2); %ntpts
D=size(t,1); %numsensors

if(nargin<2)
  Q=D-1;
end;

%%%%%
%% initialise using pca
[pcco,m_w]=pca(t',Q);

m_w_st=m_w;

% d sensors
% n repeats
% t=Wx+mu+e
% W  dxq
% mu dx1
% t  dxn
% x  qxn
% Cov(t)=WW'+sigma_e^2 I

m_x=pinv(m_w)*t;

m_x=randn(size(m_x));

tau=1./var(squash(t-m_w*m_x));

aa_0=0.001;ba_0=0.001;
  
at_0=0.001;bt_0=0.001;
  
alphamax=(aa_0+D/2)/ba_0;
alphamin=(aa_0+D/2)/(ba_0+(var(squash(m_w_st(:,1)))*pcco(1))^2*D/2);

% initial value for alpha:
%as=sqrt(alphamax*alphamin);
as=alphamin;

alpha=ones(Q,1)*as;

Sigma_w=eye(Q);

%C=m_w*m_w'+(1./tau).*eye(D);

WtW=Sigma_w+m_w'*m_w;

no_alpha_update=0;

eyeD=eye(D);
eyeQ=eye(Q);
zerosQ=zeros(Q);

ss_t=sum(t.^2,1);

%% iterate VB updates (from Bishop PCA paper)
ft_progress('init','gui');
for it=1:niters,
  ft_progress(it/niters);
  Sigma_x=pinv(eyeQ+tau*(WtW));
  
  m_x=(tau*Sigma_x*m_w'*t);  
  
  sumexpxxt=zerosQ;
  for n=1:N,
    %m_x(n,:)=tau*Sigma_x*m_w'*(t(:,n));  
    %expxxt(n,:,:)=Sigma_x + m_x(:,n)*m_x(:,n)';
    expxxt=Sigma_x + m_x(:,n)*m_x(:,n)';
    sumexpxxt=sumexpxxt+expxxt;
  end;
    
  %Sigma_w=pinv(diag(alpha)+tau*squeeze(sum(expxxt,1)));
  Sigma_w=pinv(diag(alpha)+tau*squeeze(sumexpxxt));
  
  %for d=1:D,
  %sumn=zeros(Q,1);
  %for n=1:N,
  % sumn=sumn+m_x(n,:)'*t(d,n);
  %end;
  % m_w(d,:)=tau*Sigma_w*sumn; 
  %m_w(d,:)=tau*Sigma_w*m_x'*t(d,:)';
  %end;
  
  m_w=(tau*Sigma_w*m_x*t')';

  WtW=Sigma_w+m_w'*m_w;

  if(~no_alpha_update)
    a_alpha=aa_0+D/2;
    for q=1:Q,
      b_alpha=ba_0+sum(m_w(:,q).^2)/2;
      alpha(q)=a_alpha/b_alpha;
    end;
  end;
  
  a_tau=at_0+(N-dofs_lost)*D/2;
  
  sumn=0;
  for n=1:N,
    expxxt=Sigma_x + m_x(:,n)*m_x(:,n)';        
    %sumn=sumn+trace(WtW*squeeze(expxxt(n,:,:)));
    sumn=sumn+trace(WtW*squeeze(expxxt));
    %sumn=sumn+sum(t(:,n).^2)+trace(WtW*squeeze(expxxt(n,:,:)))-2*t(:,n)'*m_w*m_x(n,:)';
  end;
  
  sumn=sumn+sum(ss_t-2*diag(t'*m_w*m_x)');
  
  b_tau=bt_0+sumn/2;

  tau=a_tau/b_tau;  
  
  taus(it)=tau;
end;

ft_progress('close');
  
%% compute covariance matrix
C=(1./tau).*eyeD;

for q=1:Q,
    tmp=m_x(q,:)*m_x(q,:)'/(N-dofs_lost)+Sigma_x(q,q);
    C=C+(m_w(:,q)*m_w(:,q)'+Sigma_w(q,q))*tmp;
end;

%% renormalise output
C=nor^2*C;
m_w=nor*m_w;
tau=1./(nor^2*1./tau);
taus=1./(nor^2*1./taus);
Sigma_w=nor^2*Sigma_w;
alpha=1./(nor^2*1./alpha);

alphamin=1./(nor^2*1./alphamin);
alphamax=1./(nor^2*1./alphamax);

