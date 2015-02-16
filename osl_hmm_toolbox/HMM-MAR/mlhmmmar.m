function [hmm,pred,gcovm] = mlhmmmar (X,T,hmm,Gamma)
% Given the state time courses estimation, does a last estimation of each MAR using (local) ML
%
% INPUT
% X             observations
% T             length of series
% hmm           HMM-MAR structure
% Gamma         p(state given X) - has to be fully defined

%
% OUTPUT
% hmm           HMM-MAR structure with the coefficients and covariance matrices updated to
%                   follow a maximum-likelihood estimation
% pred          predicted response
% gcovm         covariance matrix of the error for the entire model
%
% Author: Diego Vidaurre, OHBA, University of Oxford

order = hmm.train.order;
timelag = hmm.train.timelag;
exptimelag = hmm.train.exptimelag;
orderoffset = hmm.train.orderoffset;
zeromean = hmm.train.zeromean;
[orders,order] = formorders(order,orderoffset,timelag,exptimelag);
Sind = formindexes(orders,hmm.train.S)==1; S = hmm.train.S==1; 
if ~hmm.train.zeromean, Sind = [true(1,size(X,2)); Sind]; end
regressed = sum(S,1)>0;

K = size(Gamma,2);
ndim = size(X,2);

R =  getresiduals(X,T,Sind,order,orderoffset,timelag,exptimelag,zeromean);
XX = formautoregr(X,T,orders,order,hmm.train.zeromean);
Sind = Sind==1;

pred = zeros(size(R));
for k=1:K
    if isfield(hmm.state(k).W,'S_W'), hmm.state(k).W = rmfield(hmm.state(k).W,'S_W'); end
    if all(hmm.train.S==1)
        hmm.state(k).W.Mu_W = pinv(XX .* repmat(sqrt(Gamma(:,k)),1,size(XX,2))) * R;
    else
        hmm.state(k).W.Mu_W = zeros(size(hmm.state(k).W.Mu_W));
        for n=1:ndim
            if ~regressed(n), continue; end
            hmm.state(k).W.Mu_W(Sind(:,n),n) = pinv(XX(:,Sind(:,n)) .* repmat(sqrt(Gamma(:,k)),1,sum(Sind(:,n)))) * R(:,n);
        end
    end
    predk = XX * hmm.state(k).W.Mu_W;
    pred = pred + repmat(Gamma(:,k),1,ndim) .* predk;
    e = R(:,regressed) - predk(:,regressed);
    if strcmp(hmm.train.covtype,'diag')
        hmm.state(k).Omega.Gam_rate(regressed) = 0.5* sum( repmat(Gamma(:,k),1,sum(regressed)) .* e.^2 );
    elseif strcmp(hmm.train.covtype,'full')
        hmm.state(k).Omega.Gam_rate(regressed,regressed) =  (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
        hmm.state(k).Omega.Gam_irate(regressed,regressed) = inv(hmm.state(k).Omega.Gam_rate(regressed,regressed));
    end
end
ge = R(:,regressed) - pred(:,regressed);
if strcmp(hmm.train.covtype,'uniquediag')
    hmm.Omega.Gam_rate(regressed) = 0.5* sum( ge.^2 );
elseif strcmp(hmm.train.covtype,'uniquefull')
    hmm.Omega.Gam_rate(regressed,regressed) =  (ge' * ge);
    hmm.Omega.Gam_irate(regressed,regressed) = inv(hmm.Omega.Gam_rate(regressed,regressed));
end
gcovm = (ge' * ge) / size(R,1);

