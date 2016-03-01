% OHBA_BAYES_GLM_TUTORIAL
%
% The OHBA GLM Bayes Tutorial. 
%
% Written by Giles Colclough, Cameron Higgins, Sam Harrison and Andrew
% Quinn in Februrary 2016. 
%


%	Copyright 2016 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: GilesColclough $
%	$Revision:$
%	$LastChangedDate$
%	Contact: giles.colclough@ohba.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 25-02-2016




%% Introduction

% This script is designed to give you an intro to Bayesian methods. It used a
% toy problem where classical statistics does a fair job, but the problem is
% simple, and familiar, to ease learning. 
% 
% Using the example of  GLM we go through the classical solution, a grid
% estimation technique, several sampling methods and VB. 
% 
% An EP solution is left as a reader exercise for the adventurous, and 
% mathematically gifted. (If you get code working, let me know! - I'd like it.)
%
% Work through the script a cell at a time. Do read the text, and think
% about what's going on. There may be some blanks to fill in, between the
% dotted lines. Think about the solution before you look it up / ask.




% Here are some quick helper functions. 

% prettification
tidyAxes = @(h) set(h,...
                   'FontName', 'Helvetica', ...
				   'FontSize', 16, ...
				   'Box', 'on', ...
				   'YGrid', 'on', ...
				   'XGrid', 'off', ...
				   'TickDir', 'in', ...
				   'TickLength', [0.005 0.005], ...
				   'XMinorTick', 'off', ...
				   'YMinorTick', 'off', ...
				   'XColor', [0.3 0.3 0.3], ...
				   'YColor', [0.3 0.3 0.3], ...
				   'LineWidth', 2);
			   
% computation without numerical over/underflow
log_sum_exp = @(b) max(b) + log(sum(exp(b - max(b))));

%% Let's make some fake data - or load it from a .mat file
% Some variable is related to age and age^2
nSubjects = 10;
age       = [24 25 26 29 32 33 37 39 43 48]';

% come back later and try with more subjects! 
%--------------------------------------------------------------------------
% nSubjects = 100;
% age       = sort(20 + 30*rand(nSubjects,1));
%--------------------------------------------------------------------------

% create the ground truth. y = b0 + b1 * age + b2 * age^2 + noise
xscale     = 47;
noiseScale = 0.12;
beauty     = -age.^2 ./ xscale^2 +  2 * age ./ xscale + noiseScale * trnd(10, nSubjects, 1);

%% First things first: let's plot the data!
% plot beaty against age. 
figure('Name', 'raw data', 'Color', 'w');
plot(age, beauty, 'b+', 'MarkerSize', 6);
tidyAxes(gca);
xlabel('age');
ylabel('beauty');
ylim([0.3 1.5]);

%% Quick fit a line. 
% based on the plot, you think you might try a linear or quadratic
% regression. Build design matrices for each case.
% Don't forget the intercept!
% -------------------------------------------------------------------------
X_quadratic = [ones(nSubjects,1), age, age.^2];
X_linear    = [ones(nSubjects,1), age];
% -------------------------------------------------------------------------

%% Point estimate
% what's the normal way we solve a glm?
% Find the solution to Y = X * b + e
% 
% Hint: the answer is at https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)#The_general_problem
% If you're not sure, just ask!

% -------------------------------------------------------------------------
% a number of solutions are possible in Matlab. All, in some way, compute
% the pseudoinverse.
% beta_LS = designMat \ beauty;
% beta_LS = pinv(designMat) * beauty;
% beta_LS = (designMat' * designMat) \ (designMat' * beauty);
% beta_LS = pinv(designMat' * designMat) * designMat' * beauty;

% -------------------------------------------------------------------------

% Plot your solution over the data
yhat_linear    = X_linear    * beta_linear;
yhat_quadratic = X_quadratic * beta_quadratic;

figure('Name', 'Point estimate solution', 'Color', 'w');
plot(age, beauty, 'b+', 'MarkerSize', 6);
hold on;
plot(age, yhat_linear,    'k-', 'linewidth', 2);
plot(age, yhat_quadratic, 'r-', 'linewidth', 2);
tidyAxes(gca);
xlabel('age');
ylabel('beauty');
ylim([0.3 1.5]);

% we could, if we wanted, also estimate the standard deviation of the error
% term. 
residuals = beauty - yhat_linear; 
sigma_OLS = sqrt(residuals' * residuals ./ (nSubjects - size(X_linear,2)));

%% Likelihood function
% We're now going to move into a range of Bayesian solutions. 
% It's important that you can write a function which returns the likelihood
% of the data given any set of parameters. 

% Our model is Y = X * b + e, e ~ N(0,sigma^2).
% What's the likelihood of data Y in terms of regressors X, and the
% parameters, b and sigma?

% Write an anonymous function to compute this likelihood for any b and
% sigma.
%
% Hint: it's a normally distributed! Think about the mean and variance of Y.
%
% Hint 2: if you have the likelihood for one measurement, what's the joint 
%         probability of all measurements?
%
% Tip: an anonymous function looks like this: f = @(x) display(x);
%                      and is used like this: f(3);


% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% let's plot how that likelihood looks for the linear component of our
% model,
betaTest = -0.01:0.001:0.03;
for iTest = 1:length(betaTest),
	l(iTest) = likelihood(X_linear, [beta_linear(1); betaTest(iTest)], sigma_OLS);
end
normConstant = sum(l);
figure('Name', 'example likelihood for linear component of regression', 'Color', 'w');
plot(betaTest, l ./ normConstant, 'linewidth', 2);
tidyAxes(gca);
xlabel('linear coefficient');
ylabel('probability density');

% add on the value which OLS gave us
hold on;
plot(beta_linear(2), likelihood(X_linear, beta_linear, sigma_OLS)./normConstant, ...
	 'r+', 'MarkerSize', 6);
 
% It's going to be quite helpful to have the logarithm of this likelihood
% function. Here it is!
log_likelihood = @(X, b, sigma) sum(- (beauty - X*b).^2 ./ sigma ...
	                                - 0.5 .* log(2*pi*sigma.^2));
								
%% Prior function
% We're doing bayes, so we have to specify our prior belief in the
% parameters. As a reminder, our parameters are each component of the beta
% vector in the regression, and the standard deviation of the noise. We're
% trying to estimate all of these, so we'll need a prior on all of them. 

% If we're expressing absolutely no prior belief in our parameters, a
% standard choice for beta and sigma in linear regression is 
% (beta, sigma) ~ 1/sigma. The parametrisation w.r.t sigma reflects its
% uncertainty on a logarithmic scale. If you want to know more, ask.

% write an anonymous function to return the prior probability of beta and
% sigma. Try writing the logarithm, too.

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

%% Posterior function
% write a function which returns the log posterior probability, to within
% a constant, given the log likelihood and log prior

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

%% Grid solution
% for simple cases, we can just compute the value of the posterior
% probability over a range of parameter values. We're going to visualise
% the posterior probability of our linear model using the OLS estimate of
% sigma.

% set up the parameter space
xx = -0.8:0.05:0.8;
yy = -0.03:0.002:0.05;
[b1m, b2m] = meshgrid(xx, yy);

% compute the log posterior for our linear model, using the OLS sigma
% estimate, for each value of b1m and b2m.
for i = length(yy):-1:1,
	for j = length(xx):-1:1,
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
	end%for
end%for

% visualise!
normConstant = log_sum_exp(post_plot(:));
figure('Name', 'Grid posterior solution', 'Color', 'w');
imagesc(xx, yy, post_plot - normConstant);
tidyAxes(gca);
set(gca, 'YGrid', 'off');
xlabel('Intercept');
ylabel('Linear coefficient');
colorbar

% find the maximum a posteriori solution (the square with maximum posterior
% probability)
max_p = max(post_plot(:));
[best_b2, best_b1] = find(post_plot == max_p);
hold on
plot(xx(best_b1), yy(best_b2), 'ro', 'MarkerSize', 6);

% add on the OLS solution
plot(beta_linear(1), beta_linear(2), 'k+', 'MarkerSize', 6);

% are the results close by? what would happen if we changed the priors?

% we can also plot the 95% confidence bound.
alpha                = 0.05;
[pdf_sorted,reorder] = sort(exp(post_plot(:) - normConstant));
cdf                  = cumsum(pdf_sorted);
cdf_plot             = zeros(size(post_plot));
cdf_plot(reorder)    = cdf;
contour(xx, yy, cdf_plot, [alpha, alpha], 'r--');

% that worked quite well for two parameters, when we new the rough range of
% the posterior. What if you wanted to take the same approach for our
% quadratic model, with three parameters in the beta vector, and wanted to
% include inference about sigma? How would you go about displaying the
% results? How much memory would you need? What if you had more parameters
% still?

%% Samplers!
% Now we'll try out some sampling methods. If you remember, MCMC sampling
% methods generate a chain of dependent, locally correlated random points
% in the parameter space. After a certain amount of time, they are
% guaranteed to converge to the target (posterior) distribution. How long
% will they take? Good question...

%% 1. MH sampler
% we'll use the metropolis-hastings sampler provided in matlab. There are
% other samplers available in netlab, a matlab package, and it's not hard
% to write your own. Most major software packages (R, python, julia) will
% have samplers available.

% we need to slightly alter our posterior function, to accept a single
% vector of parameters. We'll put log(sigma) as the first entry, and the
% b-vector next.
log_post_sample = @(theta) log_posterior(X_quadratic, theta(2:end)', exp(theta(1)));

% we need to provide some starting values, and set some parameters for the
% sampler. 
theta_init = [log(0.1),   0,   0.05, -0.001; ... % chain 1
	          log(0.05),  0.1  0,     0.001; ... % chain 2
			  log(0.15), -0.1, 0,     0];        % chain 3
nSamples   = 50000;
nWarmup    = 0; % or, burn-in. Set to zero, here, for illustration.
nChains    = 3;

% we need to provide a proposal distribution. This describes how we jump
% from one sample to the next. For ease of use, we'll just use a normal
% distribution, but we'll allow the scale to be different for different
% parameters. Play around with these to see how it affects the inference. 
%
% NB: we'll move aroun in the logarithmic space of sigma, to ensure that
% sigma is always inferred as a positive variable
%
% In more sophisticated samplers, these scalings can be tuned from early
% runs of the chain, and can incorporate correlations between variables. 
%
% how did I do it? I picked a few values, reviewed, picked a few more. I
% then based the scale on the covariance of the posterior I was estimating.
% I'm not sure I did a very good job!
scales      = [0.3 0.05 0.0005 0.00002]; % or try: [0.3 0.05 0.005 0.0002]; 
covariance  = diag(scales.^2);

% you can come back later and try this setup, once you've run it once.
%--------------------------------------------------------------------------
% scales      = [0.3 0.1 0.008 1e-4];
% covariance  = diag(scales) * [1  0     0     0;     ...
% 	                          0  1    -0.99  0.98;  ...
% 							  0 -0.99  1    -0.995; ...
% 							  0  0.98 -0.995 1]     * diag(scales);
%--------------------------------------------------------------------------

% a couple of quick intermediary computations
R           = chol(covariance);
log_det_cov = 2 * sum(log(diag(R)));

% We need to generate a new sample, when we're at position x. 
% matlab doesn't include a multivariate normal sampler, so there's a bit of
% trickery coming up.
p = size(R,1);
sample_generator = @(x) x + (R' * randn(p,1))';

% we also need to define the log-pdf of a move to position y, from position x.
% this is just the multivariate normal pdf with mean x.
% log_q (y,x) = - 0.5 log(det(covariance)) exp(-0.5 (y - x) * inv(covariance) * (y - x)');
log_prop_pdf = @(y,x) - 0.5 .* sum((R' \ (y - x)').^2) - 0.5 * log_det_cov;

% let's run the sampler
clear theta_sample
parfor iChain = 1:nChains,
[theta_sample(:,:,iChain), acceptRate(iChain)] = mhsample(theta_init(iChain, :), nSamples, ...
	                                                      'logpdf',     log_post_sample,   ...
														  'logproppdf', log_prop_pdf,      ...
														  'proprnd',    sample_generator,  ...
														  'burnin',     nWarmup,           ...
														  'nchains',    1);
end%parfor

% and plot the output!
figure('Name', 'MH sampler chains');
variableNames = {'log \sigma', 'b_1', 'b_2', 'b_3'};
for iVar = 1:4,
subplot(4,1,iVar);
plot(squeeze(theta_sample(:,iVar,:)));
tidyAxes(gca);
ylabel(variableNames{iVar});
end%for

figure('Name', 'MH marginal posterior distributions', 'Color', 'w');
nBins = min([ceil(nSamples./10), 200]);
for iVar = 1:4,
subplot(1,4,iVar);
tmp = theta_sample(:,iVar,:);
[n, x] = hist(tmp(:), nBins);
width  = 0.8;
hh     = bar(gca, x, n, width);
set(hh, 'FaceColor', [120, 120, 120]/255, ... % set to be grey
		'EdgeColor', 'none');
tidyAxes(gca);
xlabel(variableNames{iVar});
set(gca,...
	'YGrid', 'off', ...
	'box',   'off', ...
	'YColor', get(gca, 'Color'));
end%for

% Notice several things about the result. 
% 1. It's awful! Why? Have a look at the correlation of the samples you've
% drawn:
display('Correlations between variables: ');
display(corr(theta_sample(:,:,1)));

% 2. There is strong bias from the iniial samples. Until the sampler has
% reached the target distribution, the initial path strongly biases the
% results. For this reason, it is usual to throw away the first 50% of
% samples as a 'warmup' period. The visual check on the chains is
% important, too. 

% 3. Chains show really poor mixing. It's clear they're not exploring the
% same space. It's important to start your chains from different positions,
% to help identify this problem. 

% The massive correlations between parameters make it really difficult for 
% the sampler to properly explore the space. What's the solution? Build in
% correlations to the MH jumping kernel, by changing the covariance we
% specified above. I've provided an example above you can try out.


%%% The section below can be left out.
%{
%% 2. Gibbs sampler
% The Gibbs sampler draws one parameter at a time. While correlations make 
% it slow to converge, it can cope (eventually) with such posteriors.
% we'll set up an mh sampler that draws each parameter independently. 
%
% Often, Gibbs samplers are used when you can derive by hand the
% appropriate distribution for each individual parameter. You can, however,
% just use the metropolis sampler again, but moving about on each parameter
% in turn. That's what we'll try here.
%
% The code will look a bit more complicated, but in essence it's simpler!
% We define 1D normal distributions on which to jump in each parameter. And
% we take it in turns to update each parameter, conditional on all of the
% others.
% we'll need normal pdf functions for a single variable
log_norm_pdf = @(y,x,scale) - 0.5 * (y - x)^2 ./ scale^2 - 0.5 * sqrt(2 * pi * scale);
norm_rand    = @(x, scale) x + scale * randn(1);
% we might want to change our scales
theta_init = [log(0.1),   0,   0.05, -0.001; ... % chain 1
	          log(0.05),  0.1  0,     0.001; ... % chain 2
			  log(0.15), -0.1, 0,     0];        % chain 3
scale    = [0.3 0.4 0.01 0.0007];  
nSamples = 5000;
theta_Gibbs    = zeros(nSamples, 4, nChains);
acceptanceRate = zeros(4,1);
parfor iChain = 1:nChains,
	% set up first sample
	theta = theta_init(iChain,:);
	for iSample = 2:nSamples,
		% a bit of progress reporting
		if ~mod(iSample, 1000),
			fprintf('Gibbs sample %d from %d in chain %d. \n', iSample, nSamples, iChain);
		end%if
		
		% note how the value from the last step in the chain gets passed in.
		
		% s.d.
		theta(1) = mhsample(theta(1), 1,                                                   ...
	                     'logpdf',     @(s) log_posterior(X_quadratic, theta(2:end)', exp(s)),  ...
						 'logproppdf', @(y,x) log_norm_pdf(y,x,scale(1)),                  ...
						 'proprnd',    @(x) norm_rand(x,scale(1)),                         ...
						 'burnin',     0,                                                  ...
						 'nchains',    1);
		% intercept
		theta(2) = mhsample(theta(2), 1,                                                   ...
	                     'logpdf',     @(b) log_posterior(X_quadratic, [b theta(3:4)]', exp(theta(1))),  ...
						 'logproppdf', @(y,x) log_norm_pdf(y,x,scale(2)),                  ...
						 'proprnd',    @(x) norm_rand(x,scale(2)),                         ...
						 'burnin',     0,                                                  ...
						 'nchains',    1);
		% linear term
		theta(3) = mhsample(theta(3), 1,                                                   ...
	                     'logpdf',     @(b) log_posterior(X_quadratic, [theta(2) b theta(4)]', exp(theta(1))),  ...
						 'logproppdf', @(y,x) log_norm_pdf(y,x,scale(3)),                  ...
						 'proprnd',    @(x) norm_rand(x,scale(3)),                         ...
						 'burnin',     0,                                                  ...
						 'nchains',    1);
		% quadratic term
		theta(4) = mhsample(theta(4), 1,                                                   ...
	                     'logpdf',     @(b) log_posterior(X_quadratic, [theta(2:3) b]', exp(theta(1))),  ...
						 'logproppdf', @(y,x) log_norm_pdf(y,x,scale(4)),                  ...
						 'proprnd',    @(x) norm_rand(x,scale(4)),                         ...
						 'burnin',     0,                                                  ...
						 'nchains',    1);
		% collect up
		theta_Gibbs(iSample,:,iChain) = theta;
	end%for
end%parfor
% plot the output!
figure('Name', 'MH within Gibbs sampler chains');
for iVar = 1:4,
subplot(4,1,iVar);
plot(squeeze(theta_Gibbs(:,iVar,:)));
tidyAxes(gca);
ylabel(variableNames{iVar});
end
% not much better, perhaps. Lots slower, too.
%}

%% 3. Block Gibbs sampler. 
% The problem with this example is that variables in the regression are
% highly correlated. It turns out that with a bit of maths, you can
% factorise the posterior into the part with sigma, and the part with the
% regression coefficients. Sampling each part in turn, where the regression
% coefficients are drawn together, gets over this problem. 
%

% Here's a blue-peter version.
% We're not even going to run multiple chains. Perhaps you could try
% chainging the code and having a look at its convergence properties over 
% three chains?
[beta_block_Gibbs, sigma_block_Gibbs] = robust_Bayes_glm(beauty, X_quadratic, inf, nWarmup, nSamples, false);

% have a look at the chain
figure('Name', 'Block Gibbs sampler', 'Color', 'w');
iVar = 1;
subplot(4,1,iVar);
plot(log(squeeze(sigma_block_Gibbs(iVar,:)')));
tidyAxes(gca);
ylabel(variableNames{iVar});
for iVar = 2:4,
subplot(4,1,iVar);
plot(squeeze(beta_block_Gibbs(iVar-1,:)'));
tidyAxes(gca);
ylabel(variableNames{iVar});
end

% histogram the results
figure('Name', 'Estimated beta', 'Color', 'w');
nBins = min([ceil(nSamples./10), 200]);
for iPlot = 1:3,
	subplot(1,3,iPlot);
	[n, x] = hist(beta_block_Gibbs(iPlot, :), nBins);
	width  = 0.8;
    hh     = bar(gca, x, n, width);
	set(hh, 'FaceColor', [120, 120, 120]/255, ... % set to be grey
            'EdgeColor', 'none');
	tidyAxes(gca);
	set(gca,...
		'YGrid', 'off', ...
		'box',   'off', ...
		'YColor', get(gca, 'Color'));
end%for

figure('Name', 'Estimated sigma', 'Color', 'w');
[n, x] = hist(sigma_block_Gibbs, nBins);
width  = 0.8;
hh     = bar(gca, x, n, width);
set(hh, 'FaceColor', [120, 120, 120]/255, ... % set to be grey
		'EdgeColor', 'none');
tidyAxes(gca);
set(gca,...
	'YGrid', 'off', ...
	'box',   'off', ...
	'YColor', get(gca, 'Color'));

% plot some samples from the posterior onto the data
nSims = 30;
predY = X_quadratic * beta_block_Gibbs(:,randperm(nSamples,nSims));
figure('Name', 'samples from posterior', 'Color', 'w');
hold on
for iSim = 1:nSims,
	plot(age, predY(:,iSim), 'Marker', 'none', ...
		 'Color', [255 62 150]/255, 'LineWidth', 1.5);
end%for
plot(age, beauty,  'o', 'Color', 'k');
tidyAxes(gca);
xlabel('age');
ylabel('beauty');

% samples must be treated as an entire vector of parameters. You can't mix
% up the order of the linear coefficient parameters against the intercept
% parameters, for example. To see why, scatter plot samples from the
% posterior, first maintaing the order, and then scrambling one w.r.t the
% other.


% Do you have confidence in any of these parameters being non-zero? Print
% out the 95% highest-density intervals, or appropriate percentiles, for
% the relevant parameters
%--------------------------------------------------------------------------
fprintf('\n\nSigma: \n\t5%%,  \t 25%%,      50%%,      75%%,      95%%.');
display(prctile(sigma_block_Gibbs, [5 25 50 75 95]));
fprintf('Regression parameters (intercept, linear, quadratic): \n\t5%%,  \t 25%%,      50%%,      75%%,      95%%.');
display(prctile(beta_block_Gibbs', [5 25 50 75 95])');
%--------------------------------------------------------------------------

%% VB

% Now let's try to solve the same problem using Variational Bayes. This
% method solves an approximate posterior function, where the posterior is
% separated ('factorised') into independent functions of each
% variable. Then, rather than working with full computational expressions
% of probability ditributions as we did above, we just work with a small  
% number of parameters (eg mean and variance) that fully describe each 
% distribution.

% to start with, can you write out on paper an expression for the
% factorised (approximate) posterior?

% We now turn our attention to the prior.
% For VB, we are required to set priors that are of a conjugate form to the
% likelihood function. Consequently we must have a normal distribution 
% (with parameters mu and sigma) for the regression coefficients and a gamma
% distribution (with parameters b and c) for the noise precision.
clear prior;
% Normal distribution over the regression coefficients with zero mean and
% large variance:
%-------------------------------------------------------------
prior.beta.mu0=zeros(size(beta_quadratic));
prior.beta.Sigma0=10*eye(length(beta_quadratic));
%-------------------------------------------------------------
%Gamma distribution over the noixe precision with mean of one and large
%variance:
%-------------------------------------------------------------
prior.lambda.b=100;
prior.lambda.c=0.01;
%-------------------------------------------------------------

% Let's plot these priors to see what they look like.
x1=linspace(-20,20,1000);
x2=linspace(0,10,1000);
beta1_pdf=normpdf(x1,prior.beta.mu0(1),sqrt(prior.beta.Sigma0(1,1)));
beta2_pdf=normpdf(x1,prior.beta.mu0(2),sqrt(prior.beta.Sigma0(2,2)));
beta3_pdf=normpdf(x1,prior.beta.mu0(3),sqrt(prior.beta.Sigma0(3,3)));
lambda_pdf=gampdf(x2,prior.lambda.c,prior.lambda.b);
figure()
subplot(2,2,1);plot(x1,beta1_pdf); title('Y axis intercept')
subplot(2,2,2);plot(x1,beta2_pdf); title('Age Coefficient')
subplot(2,2,3);plot(x1,beta2_pdf); title('Age^2 Coefficient')
subplot(2,2,4);plot(x2,lambda_pdf); title('Noise Precision')
xlim([0,3])
suptitle('Prior distributions')
%%

% Now, let's run the VB code. This calculates new values for the paramters 
% of the factorised posterior in terms of the remaining paramters; ie the
% update for the beta terms are calculated in terms of the lambda
% distribution, and vice versa. But as we don't know either of these
% distributions to start with, we run a number of iterations, using the
% most recent estimate of one paramter to inform the estimate of the other.

%let's initialise our estimates of these parameters to the prior:
beta.mu=prior.beta.mu0;
beta.Sigma=prior.beta.Sigma0;
lambda.b=prior.lambda.b;
lambda.c=prior.lambda.c; 

%Now run successive iterations and plot each time:
figure(); suptitle('Factorised Posterior Functions Over Successive Iterations');
for iteration=1:10
    % write a function beta_update that updates the structure beta given the 
    % following parameters (see slides for the update equations)
    %-------------------------------------------------------
    beta=beta_update(prior.beta,lambda,beauty,X_quadratic);
    %-------------------------------------------------------
    
    %plot update for beta:
    x1=linspace(beta.mu(1)-4*beta.Sigma(1,1),beta.mu(1)+4*beta.Sigma(1,1),1000);
    beta1_pdf=normpdf(x1,beta.mu(1),beta.Sigma(1,1));
    subplot(2,2,1);hold on;  plot(x1,beta1_pdf); %xlim([x1(1),x1(end)]);
    ylabel('Q(\beta_1)');xlabel('\beta_1');title('Y axis intercept')
    %Plot OLS estimate as a line:
    %line([beta_quadratic(1),beta_quadratic(1)],[0 max(beta1_pdf)],...
     %       'LineWidth',2,'Color','r');
        
    x2=linspace(beta.mu(2)-4*beta.Sigma(2,2),beta.mu(2)+4*beta.Sigma(2,2),1000);
    beta2_pdf=normpdf(x2,beta.mu(2),beta.Sigma(2,2));
    subplot(2,2,2); hold on; plot(x2,beta2_pdf); %xlim([x2(1),x2(end)]);
    ylabel('Q(\beta_2)');xlabel('\beta_2');title('Age Coefficient')
        
    x3=linspace(beta.mu(3)-4*beta.Sigma(3,3),beta.mu(3)+4*beta.Sigma(3,3),1000);
    beta3_pdf=normpdf(x3,beta.mu(3),beta.Sigma(3,3));
    subplot(2,2,3); hold on; plot(x3,beta3_pdf); %xlim([x3(1),x3(end)]);
    ylabel('Q(\beta_3)');xlabel('\beta_3');title('Age^2 Coefficient')
    
    
    % now write a function lambda_update that updates the structure LAMBDA
    % given the following parameters (see slides for the update equations)
    %-------------------------------------------------------------
    lambda = lambda_update(prior.lambda,beta,beauty,X_quadratic);
    %-------------------------------------------------------------
    
    %and plot this update:
    x3=linspace(0,lambda.b*lambda.c+lambda.b^2*lambda.c,1000);
    lambda_pdf=gampdf(x3,lambda.c,lambda.b);
    sigma_cdf=1-gamcdf(x3,lambda.c,lambda.b);
    sigma_pdf=-diff(1-gamcdf(x3,lambda.c,lambda.b)); %convert from precision back to std
    subplot(2,2,4); hold on; plot(x3(2:end).^-0.5,sigma_pdf); xlim([0,x3(2).^-0.5])
    ylabel('Q(\sigma)');xlabel('\sigma');title('Noise Std')
    
    pause(1);
end

%% 
% This has shown us the individual factorised posterior estimates, but how
% do these differ from the true posterior? We will now plot the noise
% precision against the intercept estimate - how do you expect this to
% differ from the true posterior plotted in the grid solution?

% compute the log posterior for our linear model, using the OLS sigma
% estimate, for each value of b1m and b2m.
for i = 1:length(sigma_pdf)
	for j = 1:length(beta1_pdf)
% -------------------------------------------------------------------------
		post_plot(i,j) = beta1_pdf(j)*sigma_pdf(i);
% -------------------------------------------------------------------------
	end%for
end%for

% visualise!
normConstant = log_sum_exp(post_plot(:));
figure('Name', 'Grid posterior solution', 'Color', 'w');
imagesc(x3(2:end).^-0.5, x1, post_plot);
tidyAxes(gca);
set(gca, 'YGrid', 'off');
xlabel('Noise std');
ylabel('Intercept');
colorbar


%% All done?
% Try generating more data - use 100 subjects this time. Does that improve
% the certainty of your estiamtes? Grab some weather data or climate change
% data. Can you fit a model with a linear trend, and some seasonal /
% monthly variation?

%% All too much? 
% Never fear. There are excellent packages which don't even
% require you to set a proposal distribution. So long as you can write your
% likelihood and prior in maths, (maybe also as a function, as you've done
% here), they'll work for you. They key is to be able to assess their
% performance.
%
% Top sampling packages are pymc3 (written in python) and STAN (separate,
% v. easy syntax, super clever NUTS sampler, compiles itself so super fast, and great
% unless you have tons of parameters or megatons of data).
% 
% When you take on a new sampling package, try out the demos first. Obvs.

%% VB in practice. 
% VB takes a lot of maths to get right. When it works, it's fast, and can
% be very accurate. It's easy to monitor convergance, too. 
%
% However, most of the time, you're going to want someone to write the VB
% code for you. 
%
% There are easy-to-use out of the box packages, notably "infer.net". In our
% experience, they can't handle lots of data (MEG datasets, for example),
% and run on Microsoft / anything that can handle .NET frameworks.

%% example stan code for linear regression:
% (Stan needs to be installed. It does interface with matlab, R, julia and
% python).
%{
/* glm.stan */
//-------------------------------------------------------------------------
data {
    int    <lower=1>                         nSubjects; 
	int    <lower=1>                         nPredictors;
	vector <lower=0>[nSubjects]              y;
	matrix          [nSubjects, nPredictors] X;
}
//-------------------------------------------------------------------------
parameters {
	vector [nPredictors]     beta;
	real   <lower=0,upper=5> sigma; // implicit uniform prior
}
//-------------------------------------------------------------------------
model {
    /* prior */
	beta ~ normal(0, 5);

	/* likelihood */
	y ~ normal(X * beta, sigma^2);
}
//-------------------------------------------------------------------------
%}
% and to run it:
%{
#!/bin/bash
# script to run glm.stan

# Compile stan model
echo "Compiling model"
#make -f ${STAN_HOME}/makefile -j4 glm

# Set four chains running
echo "setting chains running: "
for i in {1..4}
do
    echo "$i"
    cmd="./glm sample num_warmup=1000 num_samples=2000 data file=my_data.Rdump output file=sample${i}.csv id=${i}"
    ${cmd} &
done

# summarise output
${STAN_HOME}/bin/print sample1.csv sample2.csv sample3.csv sample4.csv
echo "Done!"
%}
