function [m_w, alpha, C, invC, taus, tau, m_x, Sigma_w, Sigma_x, Q, N, D,alphamin,alphamax, maxw]=bayes_pca_ss(x,Q,niters,dofs_lost)

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

% clear
% T=1000;N=200;alpha=2/N*randn(N,N);alpha=alpha*alpha';alpha=alpha-diag(diag(alpha));x=mvnrnd(zeros(N,1),inv(eye(N)-alpha),T)';cov(x');[m_w, alpha, C]=bayes_pca2(x,100);

if(nargin<3)
    niters=30;
    dofs_lost=1;
end;

t=demean(x,2);

Q=min(Q,size(t,2));

nor=std(squash(t));
t=t/nor;

%[pcco,m_w]=pca(t',Q);V= m_w(:,1:Q);figure;plot(V(:,1));

%[VT E]  = spm_svd(t*t',exp(-8));figure;plot(VT(:,1));

[V s]  = svd(t*t');               % get temporal modes
V    = V(:,1:Q);                      %  modes to keep

%figure;plot(V(:,1));

told=t;
V=V';
t=V*t;

N=size(t,2); %ntpts
D=size(t,1); %numsensors

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
%aa_0=0.0001;ba_0=0.0001;
  
at_0=0.001;bt_0=0.001;
%at_0=0.0001;bt_0=0.0001;
%at_0=0.0000000001;bt_0=0.0000000001;
  
alphamax=(aa_0+D/2)/ba_0;
%alphamin=(aa_0+D/2)/(ba_0+max(squash(m_w_st(:,1))*pcco(1))^2*D/2)
alphamin=(aa_0+D/2)/(ba_0+(var(squash(m_w_st(:,1)))*pcco(1))^2*D/2);

as=(alphamin+alphamax)/2
as=alphamin;

alpha=ones(Q,1)*as;

Sigma_w=eye(Q);

%C=m_w*m_w'+(1./tau).*eye(D);

WtW=Sigma_w+m_w'*m_w;

%expxxt=zeros(N,Q,Q);
no_alpha_update=0;

eyeD=eye(D);
eyeQ=eye(Q);
zerosQ=zeros(Q);

ss_t=sum(t.^2,1);

tauold=tau;
  
progress('init','etf','Computing Bayes PCA');
for it=1:niters,
   
  progress(it/niters);
  perctaudiff=100*abs(tau-tauold)/tauold;
  pca_dim = min(find(alpha>((alpha(1)+alphamax)/2)))-1; 

  disp(['%diff tau=' num2str(perctaudiff) ', pca_dim=' num2str(pca_dim)]);

  if(it>5 && perctaudiff<0.5), break; end;
  
  tauold=tau;
   
  Sigma_x=pinv(eyeQ+tau*(WtW));
  
  m_x=(tau*Sigma_x*m_w'*t);  
  
  sumexpxxt=zerosQ;
  for n=1:N,
   % m_x(n,:)=tau*Sigma_x*m_w'*(t(:,n));  
    %expxxt(n,:,:)=Sigma_x + m_x(:,n)*m_x(:,n)';
    expxxt=Sigma_x + m_x(:,n)*m_x(:,n)';
    sumexpxxt=sumexpxxt+expxxt;
  end;
    
%  Sigma_w=pinv(diag(alpha)+tau*squeeze(sum(expxxt,1)));
  Sigma_w=pinv(diag(alpha)+tau*sumexpxxt);
  
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
    sumn=sumn+trace(WtW*expxxt);
    %sumn=sumn+sum(t(:,n).^2)+trace(WtW*squeeze(expxxt(n,:,:)))-2*t(:,n)'*m_w*m_x(n,:)';
  end;
  
  sumn=sumn+sum(ss_t-2*diag(t'*m_w*m_x)');
  
  b_tau=bt_0+sumn/2;
%  if(it>18)keyboard;end;

  tau=a_tau/b_tau;  
    
  taus(it)=tau;
  
end;
progress('close');

%C=m_w*(m_x*m_x'+Sigma_x)*m_w'+(1./tau).*eye(D);

%C=m_w*m_w'+(1./tau).*eye(D);
Cproj=(1./tau).*eyeD;

for q=1:Q,
    tmp=m_x(q,:)*m_x(q,:)'/(N-dofs_lost)+Sigma_x(q,q);
    Cproj=Cproj+(m_w(:,q)*m_w(:,q)'+Sigma_w(q,q))*tmp;
    %Cproj=Cproj+(m_w(:,q)*m_w(:,q)');
end;

% put back into full space
pinvV=pinv(V);

C=pinvV*Cproj*pinvV';

maxw = min(find(alpha>((alpha(1)*alphamax)^(1/2))))-1; 

invC=V(1:maxw,:)'*pinv(Cproj(1:maxw,1:maxw))*V(1:maxw,:);

Sigma_w=pinvV*Sigma_w*pinvV';
m_w=pinvV*m_w;

% renorm

C=nor^2*C;
invC=invC/nor^2;

m_w=nor*m_w;
tau=1./(nor^2*1./tau);
taus=1./(nor^2*1./taus);
Sigma_w=nor^2*Sigma_w;
alpha=1./(nor^2*1./alpha);

alphamin=1./(nor^2*1./alphamin);
alphamax=1./(nor^2*1./alphamax);

%cov(t')
