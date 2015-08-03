function [Sig_save,C_save,lambda_save] = BayesGLasso_Columnwise(S,n,Sig,C,a_lambda,b_lambda,burnin,nmc)
% Efficient Bayesian Graphical Lasso MCMC sampler using data-augmented
% block (column-wise) Gibbs sampler
%Input:
%     S = Y'*Y : sample covariance matrix * n
%     n: sample size
%     lambda:   |C|^(n/2) exp(-1/2 (SC) -  lambda/2 ||C||_1 );
%     Sig,C: initial Covariance and precision matrix C = inv(Sig);
%     burnin, nmc : number of MCMC burnins and saved samples
%     lambda ~ Ga(a_lambda,b_lambda) - shape a_lambda, scale 1/b_lambda

%Output:
%     Sig_save: p by p by nmc matrices of saved posterior samples of covariance
%     C_save: p by p by nmc matrices of saved posterior samples of precision
%     lambda: 1 by nmc vector of saved posterior samples of lambda

%  Ref:  Wang 2012 Bayesain Graphical Lasso and Efficient Posterior
%  Computation

%  Written by Hao Wang & U of South Carolina

[p] = size(S,1); indmx = reshape([1:p^2],p,p); 
upperind = indmx(triu(indmx,1)>0); 

indmx_t = indmx';
lowerind = indmx_t(triu(indmx_t,1)>0); 


C_save = zeros(p,p,nmc); Sig_save = C_save;
lambda_save = zeros(1,nmc);
 tau = zeros(p);

ind_noi_all = zeros(p-1,p);
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       
       ind_noi_all(:,i) = ind_noi;
end

    apost = a_lambda + p*(p+1)/2; 

ft_progress('init', 'text', '');
cleanProgressDisplay = onCleanup(@() ft_progress('close'));
for iter = 1: burnin+nmc    
            
       if(mod(iter,1000)==0),
           ft_progress(iter / (burnin+nmc),              ...
                       '    MCMC iter = %d out of %d\n', ...
                       iter, burnin+nmc);
       end%if

       
    % %%% Sample lambda 
    bpost  = b_lambda + sum(abs(C(:)))/2;    
    lambda = ROInets.randgamma(apost, 1.0/bpost); % replaced gamrnd(apost,1/bpost,1); using Lightspeed Mex file. Uses properties of scale 1/bpost: see code in randgamma, gamrnd, randg, and http://en.wikipedia.org/wiki/Gamma_distribution#Scaling
    
%%% sample tau off-diagonal        
    Cadjust = max(abs(C(upperind)),10^-6);        
    lambda_prime = lambda^2;  
    mu_prime = min(lambda./Cadjust,10^12);
                

    tau_temp =  1./rand_ig(mu_prime,lambda_prime);
    tau(upperind) = tau_temp;
    tau(lowerind) = tau_temp;
    

%%% sample Sig and C = inv(Sig)        
    for i = 1:p
        
      ind_noi = ind_noi_all(:,i);
 
      tau_temp = tau(ind_noi,i);
       
      Sig11 = Sig(ind_noi,ind_noi); Sig12 = Sig(ind_noi,i);
      
      invC11 = Sig11 - Sig12*Sig12'/Sig(i,i);
      
      Ci = (S(i,i)+lambda)*invC11+diag(1./tau_temp);
  
        
        Ci_chol = chol(Ci);
        
        mu_i = -Ci\S(ind_noi,i);

        beta = mu_i+ Ci_chol\randn(p-1,1);
        
        C(ind_noi,i) = beta;
        C(i,ind_noi) = beta;
        gam = ((S(i,i)+lambda)\2) .* randgamma(n/2+1); % replaces gamrnd(n/2+1,(S(i,i)+lambda)\2);
        
        C(i,i) = gam+beta'*invC11*beta;
        
        %% Below updating Covariance matrix according to one-column change of precision matrix
        invC11beta = invC11*beta;
        
        Sig(ind_noi,ind_noi) = invC11+invC11beta*invC11beta'/gam;
        Sig12 = -invC11beta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;
    end

   if iter >burnin           
        Sig_save(:,:,iter-burnin) = Sig; 
        C_save(:,:,iter-burnin) = C;
        lambda_save(iter-burnin) = lambda;
   end%if
end
end%BayesGlassoColumnwise




function y = rand_ig(theta,chi)
%RAND_IG
% THE INVERSE GAUSSIAN DISTRIBUTION
%
% The Inverse Gaussian distribution is left skewed distribution whose
% location is set by the mean with the profile determined by the
% scale factor.  The random variable can take a value between zero and
% infinity.  The skewness increases rapidly with decreasing values of
% the scale parameter.
%
%
% pdf(y) = sqrt(chi/(2*pi*y^3)) * exp(-chi./(2*y).*(y/theta-1).^2);
% cdf(y) = normcdf(sqrt(chi./y).*(y/theta-1)) + ...
%            exp(2*chi/theta)*normcdf(sqrt(chi./y).*(-y/theta-1));
%
%   where  normcdf(x) = 0.5*(1+erf(y/sqrt(2))); is the standard normal CDF
%         
% Mean     = theta;
% Variance = theta^3/chi;
% Skewness = sqrt(9*theta/chi);
% Kurtosis = 15*mean/scale;
% Mode = theta/(2*chi)*(sqrt(9*theta^2+4*chi^2)-3*theta);
%
% PARAMETERS:
%  theta - location; (theta>0)
%  chi - scale; (chi>0)
%
% SUPPORT:
%  y,  y>0
%
% CLASS:
%   Continuous skewed distribution
%
% NOTES:
%   1. There are several alternate forms for the PDF, 
%      some of which have more than two parameters
%   2. The Inverse Gaussian distribution is often called the Inverse Normal
%   3. Wald distribution is a special case of The Inverse Gaussian distribution
%      where the mean is a constant with the value one.
%   4. The Inverse Gaussian distribution is a special case of The Generalized
%        Hyperbolic Distribution
%
% USAGE:
%   randraw('ig', [theta, chi], sampleSize) - generate sampleSize number
%         of variates from the Inverse Gaussian distribution with 
%         parameters theta and chi;
%   randraw('ig') - help for the Inverse Gaussian distribution;
%
% EXAMPLES:
%  1.   y = randraw('ig', [0.1, 1], [1 1e5]);
%  2.   y = randraw('ig', [3.2, 10], 1, 1e5);
%  3.   y = randraw('ig', [100.2, 6], 1e5 );
%  4.   y = randraw('ig', [10, 10.5], [1e5 1] );
%  5.   randraw('ig');
% 
% SEE ALSO:
%   WALD distribution
% END ig HELP END inversegauss HELP  END invgauss HELP 

% Method:
%
% There is an efficient procedure that utilizes a transformation
% yielding two roots.
% If Y is Inverse Gauss random variable, then following to [1]
% we can write:
% V = chi*(Y-theta)^2/(Y*theta^2) ~ Chi-Square(1),
%
% i.e. V is distributed as a chi-square random variable with
% one degree of freedom.
% So it can be simply generated by taking a square of a
% standard normal random number.
% Solving this equation for Y yields two roots:
%
% y1 = theta + 0.5*theta/chi * ( theta*V - sqrt(4*theta*chi*V + ...
%      theta^2*V.^2) );
% and
% y2 = theta^2/y1;
%
% In [2] showed that  Y can be simulated by choosing y1 with probability
% theta/(theta+y1) and y2 with probability 1-theta/(theta+y1)
%h
% References:
% [1] Shuster, J. (1968). On the Inverse Gaussian Distribution Function,
%         Journal of the American Statistical Association 63: 1514-1516.
%
% [2] Michael, J.R., Schucany, W.R. and Haas, R.W. (1976).
%     Generating Random Variates Using Transformations with Multiple Roots,
%     The American Statistician 30: 88-90.

sampleSize = max(length(theta),length(chi));

chisq1 = randn(sampleSize,1).^2;
y = theta + 0.5.*theta./chi.*(theta.*chisq1 - ...
     sqrt(4.*theta.*chi.*chisq1 + theta.^2.*chisq1.^2) );

l = rand(sampleSize,1)>= theta./(theta+y);   

if any(l),
    y(l) = theta(l).^2./y(l);
end
end%rand_ig
% [EOF]