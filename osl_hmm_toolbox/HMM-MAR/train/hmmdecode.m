function [vpath]=hmmdecode(X,T,hmm,residuals)
%
% Viterbi and single-state decoding for hmm
% The algorithm is run for the whole data set, including those whose class
% was fixed. This means that the assignment for those can be different.
%
% INPUT
% X      Observations
% T             length of series
% hmm       hmm data structure
% residuals     in case we train on residuals, the value of those.
%
% OUTPUT
% vpath(i).q_star    maximum likelihood state sequence
%
% Author: Diego Vidaurre, OHBA, University of Oxford

tiny=exp(-700);
N = length(T);
K=hmm.K;
P=hmm.P;
Pi=hmm.Pi;

[~,order] = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);


for in=1:N
    
    q_star = ones(T(in)-order,1);
    
    alpha=zeros(T(in),K);
    beta=zeros(T(in),K);
    
    % Initialise Viterbi bits
    delta=zeros(T(in),K);
    psi=zeros(T(in),K);
    
    if in==1, t0 = 0; s0 = 0;
    else t0 = sum(T(1:in-1)); s0 = t0 - order*(in-1);
    end
    
    B = obslike(X(t0+1:t0+T(in),:),hmm,residuals(s0+1:s0+T(in)-order,:));
    
    scale=zeros(T(in),1);
    % Scaling for delta
    dscale=zeros(T(in),1);
    
    alpha(1+order,:)=Pi(:)'.*B(1+order,:);
    scale(1+order)=sum(alpha(1+order,:));
    alpha(1+order,:)=alpha(1+order,:)/(scale(1+order)+tiny);
    
    delta(1+order,:) = alpha(1+order,:);    % Eq. 32(a) Rabiner (1989)
    % Eq. 32(b) Psi already zero
    for i=2+order:T(in)
        alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
        scale(i)=sum(alpha(i,:));
        alpha(i,:)=alpha(i,:)/(scale(i)+tiny);
        
        for k=1:K,
            v=delta(i-1,:).*P(:,k)';
            mv=max(v);
            delta(i,k)=mv*B(i,k);  % Eq 33a Rabiner (1989)
            if length(find(v==mv)) > 1
                % no unique maximum - so pick one at random
                tmp1=find(v==mv);
                tmp2=rand(length(tmp1),1);
                [~,tmp4]=max(tmp2);
                psi(i,k)=tmp4;
            else
                psi(i,k)=find(v==mv);  % ARGMAX; Eq 33b Rabiner (1989)
            end
        end;
        
        % SCALING FOR DELTA ????
        dscale(i)=sum(delta(i,:));
        delta(i,:)=delta(i,:)/(dscale(i)+tiny);
    end;
    
    % Get beta values for single state decoding
    beta(T(in),:)=ones(1,K)/scale(T(in));
    for i=T(in)-1:-1:1+order
        beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i);
    end;
    
    xi=zeros(T(in)-1-order,K*K);
    for i=1+order:T(in)-1
        t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
        xi(i-order,:)=t(:)'/sum(t(:));
    end;
    
    delta=delta(1+order:T(in),:);
    psi=psi(1+order:T(in),:);
    
    % Backtracking for Viterbi decoding
    id = find(delta(T(in)-order,:)==max(delta(T(in)-order,:)));% Eq 34b Rabiner;
    q_star(T(in)-order) = id(1);
    for i=T(in)-1-order:-1:1,
        q_star(i) = psi(i+1,q_star(i+1));
    end
    
    vpath(in).q_star = q_star;
    
end

