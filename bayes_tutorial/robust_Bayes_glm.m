function [beta, sigma] = robust_Bayes_glm(Y, X, nu, nWarmup, nSamples, doPlots)
%ROBUST_BAYES_GLM  Robust Bayesian linear regression
%
% [BETA, SIGMA] = ROBUST_BAYES_GLM(Y, X, NU, NW, NS) samples of BETA and
%   SIGMA from a robust bayesian linear regression of Y = X*BETA + e, where
%   e is t_NU distributed with standard deviation SIGMA. The degrees of
%   freedom parameter controls the extent to which outliers are tolerated.
%   Set NU = 4 for standard robust regression, and NU = INF to use a normal
%   noise model. 
%
%   BETA and SIGMA contain NW (nWarmup) + NS (nSamples) samples in their
%   columns. 
%
%   X must be a design matrix for linear regression, of the same number of
%   rows (observations) as Y, and with a predictor in each column. No bias
%   term is included, if you want to model a constant, set one column of X
%   to be entirely ones. 
%
%   Y can be a single set of observations, or draws from the posterior
%   distributions of Y given some other data - e.g. if Y are samples from
%   some other fit. Each sample of Y occupies a column. 
%
%   Prior distributions are (beta, log(sigma)) ~ 1. This non-informative
%   prior should be sufficient if there are quite a few observations (say,
%   over 20 for univariate linear regression) and not too many predictors. 
%
%   Choosing NU to be Inf will yield a model for which the maximum of the
%   posterior on beta and sigma coincides with the results of standard 
%   classical linear regression. 
%
% [...] = ROBUST_BAYES_GLM(..., true) outputs example plots of regression
%   lines relative to each predictor in X, plotted against the mean of
%   input data samples Y for comparison. Also output are histograms of the
%   fitted parameters, using samples from after the warmup period only. 
%
%   This code runs one sampling chain only. It is often advised to run
%   multiple chains of sampling, and to assess convergence using the
%   Gelman-Rubin statistic. 

%	References:
%	Gelman, A. et al., "Bayesian Data Analysis, 3rd Ed.", CRC Press, 2014. 
%	

%	Copyright 2015 OHBA
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
%	$Revision: 763 $
%	$LastChangedDate: 2015-10-21 11:52:19 +0100 (Wed, 21 Oct 2015) $
%	Contact: giles.colclough@magd.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 02-Mar-2015 20:27:04

%% Check inputs
% Y can be a set of samples from marginal posterior values Y conditional on
% data
[nPoints, nSamplesIn] = size(Y);

[checkSize, kPredictors] = size(X);
assert(isequal(nPoints, checkSize),      ...
       [mfilename ':UnequalDataPoints'], ...
       'Number of data points in X and Y must be identical. \n');

% If nu is non-infinite, run robust regression
if isinf(nu),
    isRobust = false;
else
    assert(isscalar(nu) && ~mod(nu,1) && nu >= 3, ...
           [mfilename ':badNu'],                  ...
           'Nu must be a scalar degrees of freedom > 2, or infinite. \n');
    isRobust = true;
end%if

% If number of samples requested is greater than number of samples of Y,
% then sample with replacement. 
% If only one sample is provided for y (as in, measured data, not
% estimated) then flag this also
nSamplesTotal = nWarmup + nSamples;
if 1 == nSamplesIn,
    isYsampled            = false;
    sampleWithReplacement = true;
elseif nSamplesIn < nSamplesTotal,
    isYsampled            = true;
    sampleWithReplacement = true;
else
    isYsampled            = true;
    sampleWithReplacement = false;
end%if

if nargin < 6 || ~exist('doPlots', 'var'),
    doPlots = false;
end%if
    
%% Initialisation
% initialise variance parameters to be within the std dev of y.
sigmasq = var(Y(:,1));
V       = repmat(sigmasq, nPoints, 1);

% declare memory
beta  = NaN(kPredictors, nSamplesTotal);
sigma = NaN(1, nSamplesTotal);

%% Run sampling
if isRobust, 
    % must compute conditional on beta every time
    
    if ~isYsampled,
        yInd = ones(1, nSamplesTotal);
    else
        % precompute samples from y
        yInd = pick_y(nSamplesIn, nSamplesTotal, sampleWithReplacement);
    end%if
    
    for iSample = 1:nSamplesTotal,
        if ~mod(iSample, 100) || 1 == iSample,
            fprintf('%s: sample %d out of %d. \n', ...
                    mfilename, iSample, nSamplesTotal);
        end%if
        % compute conditional on beta
        Qinv      = diag(1.0 ./ V);
        VbInvChol = chol(X' * Qinv * X);
        betahat   = VbInvChol \ (VbInvChol' \ (X' * Qinv * Y(:, yInd(iSample)))); % Vb * X' * Qinv * Y(:, yInd(iSample));
        
        % use nested normal parameterisation of t-distribution
        beta(:,iSample) = betahat + VbInvChol \ randn(length(betahat), 1);
        residuals       = Y(:, yInd(iSample)) - X * beta(:,iSample);
        sSq             = residuals'*residuals ./ (nPoints - kPredictors);
        sigmasq         = rand_inv_chi_sq(nPoints - kPredictors, sSq);
        sigma(iSample)  = sqrt(sigmasq);
        V               = rand_inv_chi_sq(nu + 1,                       ...
                                          (nu * sigmasq + residuals.^2) ...
                                           ./ (nu + 1));
    end%for
    
else % normal distribution
    % precompute conditional on beta
    XTXchol   = chol(X' * X);
    
    if ~isYsampled,
        % precompute betahat
        y       = Y;
        betahat = XTXchol \ (XTXchol' \ (X' * Y)); %invXTXX * Y;
    else
        % precompute samples from y
        yInd = pick_y(nSamplesIn, nSamplesTotal, sampleWithReplacement);
    end%if
    
    for iSample = 1:nSamplesTotal,
        if ~mod(iSample, 100) || 1 == iSample,
            fprintf('%s: sample %d out of %d. \n', ...
                    mfilename, iSample, nSamplesTotal);
        end%if
        
        if isYsampled,
            y       = Y(:, yInd(iSample));
            betahat = XTXchol \ (XTXchol' \ (X' * y));
        end%if
        
        sigma(iSample)  = sqrt(sigmasq);
        beta(:,iSample) = betahat + sigma(iSample).* (XTXchol \ randn(length(betahat),1));
        
        residuals       = y - X * beta(:,iSample);
        sSq             = residuals'*residuals ./ (nPoints - kPredictors);
        sigmasq         = rand_inv_chi_sq(nPoints - kPredictors, sSq);
    end%for
end%if



%% Make some plots to show the result
if doPlots,
    % generate some simulated fits
    nSims = 30;
    predY = X * beta(:,nWarmup+randperm(nSamples,nSims));
    
    figure('Name', 'Fitted data', 'Color', 'w');
    for iPlot = 1:kPredictors,
        subplot(1,kPredictors,iPlot);
        hold on
        for iSim = 1:nSims,
            plot(X(:,iPlot), predY(:,iSim), 'Marker', 'none', ...
                 'Color', [255 62 150]/255, 'LineWidth', 1.5);
        end%for
        plot(X(:,iPlot), mean(Y,2),  'o', 'Color', 'k');
		set(gca,...
            'FontName', 'Helvetica', ...
            'FontSize', 14, ...
            'Box', 'on', ...
            'YGrid', 'off', ...
            'XGrid', 'off', ...
            'TickDir', 'in', ...
            'TickLength', [0.005 0.005], ...
            'XMinorTick', 'off', ...
            'YMinorTick', 'off', ...
            'XColor', [0.3 0.3 0.3], ...
            'YColor', [0.3 0.3 0.3], ...
            'LineWidth', 2);
    end%for
    
    % histogram the results
    figure('Name', 'Estimated beta', 'Color', 'w');
    nBins = min([ceil(nSamples./10), 200]);
    for iPlot = 1:kPredictors,
        subplot(1,kPredictors,iPlot);
        GC_histogram(beta(iPlot, (nWarmup+1):end), nBins);
        set(gca,...
            'FontName', 'Helvetica', ...
            'FontSize', 14, ...
            'Box', 'off', ...
            'YGrid', 'off', ...
            'XGrid', 'off', ...
            'TickDir', 'in', ...
            'TickLength', [0.005 0.005], ...
            'XMinorTick', 'off', ...
            'YMinorTick', 'off', ...
            'XColor', [0.3 0.3 0.3], ...
            'YColor', get(gca, 'Color'), ...
            'LineWidth', 2);
    end%for
    
    figure('Name', 'Estimated sigma', 'Color', 'w');
    GC_histogram(sigma((nWarmup+1):end), nBins);
        set(gca,...
            'FontName', 'Helvetica', ...
            'FontSize', 14, ...
            'Box', 'off', ...
            'XGrid', 'off', ...
            'YGrid', 'off', ...
            'TickDir', 'in', ...
            'TickLength', [0.005 0.005], ...
            'XMinorTick', 'off', ...
            'YMinorTick', 'off', ...
            'XColor', [0.3 0.3 0.3], ...
            'YColor', get(gca, 'Color'), ...
            'LineWidth', 2);
end%if

%% Trim out warmup
beta(:,1:nWarmup) = [];
sigma(1:nWarmup) = [];
end%robust_Bayes_glm

function ind = pick_y(nSIn, nS, sampleWithReplacement)
% pick nS samples from 1:nSIn, with or without replacement
if sampleWithReplacement,
    ind = randi(nSIn, 1, nS);
else
    ind = randperm(nSIn, nS);
end%if
end%pick_y

function sigmaSq = rand_inv_chi_sq(nu, sSq)
% sample from inverse chi-square distribution

assert(isscalar(nu),         ...
       [mfilename ':BadNu'], ...
       'Can only deal with scalar degrees of freedom. \n');
   
if ~isscalar(sSq),
    % draw random samples for each element of scale sSq
    nu = repmat(nu, size(sSq));
end%if

% sigmaSq^{-1} ~ chi_sq(nu, sSq)
% sigmaSq^{-1} ~ gamma(nu/2, nu*sSq/2)
% use wikipedia on gamma dist and help for randgamma to find:
sigmaSq = nu .* sSq ./ (2 .* randgamma(nu./2));
end%rand_inv_chi_sq
% [EOF]