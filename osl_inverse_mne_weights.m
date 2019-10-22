function [W, W_nonorm,lf] = osl_inverse_mne_weights(SensorData, LeadFields, Noise, Options)
%OSL_INVERSE_MNE_WEIGHTS creates weights for source reconstruction using MNE
%
% W = OSL_INVERSE_MNE_WEIGHTS(SensorData, LeadFields, NOISE, OPTIONS)
%   calculates cell array of weights, W, for source reconstruction. Sources
%   at dipole i are estimated with W{i}*B. 
%
% [W, W_NONORM, LF] = OSL_MNE_WEIGHTS(...) also returns the un-normalised
%   weights W_NONORM and a cell array of lead fields for each dipole, LF. 
%
% The inputs should be formatted as a set of data structures:
%
% SensorData.cov       - covariance of the sensor data (nChans x nChans)
%           .nSamples  - number of time samples used to compute covariance
%                        matrix. 
%
% LeadFields.nSources  - number of sources to estimate
%           .nDims     - dimensionality of each source (expected to be 3 at
%                        the moment)
%           .lf        - matrix of lead fields, nChans x (nDims x
%                        nSources). Fields for x, y and z components of the
%                        same dipole should appear in consecutive columns. 
%
% Noise.model          - determines specification of the noise covariance.
%                        Choices are {'empirical', 'white',
%                        'diag_datacov'}. For the first option, a noise
%                        covariance matrix must be passed in. For the
%                        second two, a scaling factor lambda should be
%                        provided, or is assumed unity. 
%      .lambda         - scaling factor on the noise covariance
%      .cov            - empirical noise covariance (nChans x nChans)
%
% Further algorithm choices are set in the subfields of the Options
% structure. 
%
% Options.ReduceRank.weights    - Boolean [true] converts to scalar weights
%                   .leadFields - Integer [3] dimensionality of lead fields
%                                 (ignored at the moment)
%
%        .sourceModel     - Choose the estimation method using a string, from:
%                           'Wens' (based on paper referenced below), 
%                           'MNE' (which uses gamma-MAP estimation of a white covariance matrix), 
%                           'MNE-INA' (integrated nested approximation for
%                           averaging over the prior on gamma)
%                           'MNE-scaled-noise' (also estimates a parameter
%                           rho for the scaling between sensor data and
%                           noise)
%                           'sparse-Bayes' (gamma-MAP estimation of a
%                           sparse Bayes solution)
%                           'rvm-beamformer' - beamformer solution using
%                           the regularised data covariance estimated using
%                           the sparse Bayes model. 
%
%        .Normalise.weights - sLoreta, norm, noise-proj, none: normalises
%                             reconstruction weights using a variety of
%                             methods. sLoreta uses normalisation from the
%                             Wens paper. 'norm' normalised the 2-norm of
%                             the weights. 'noise-proj' normalises by the
%                             projection of the noise:
%                             sqrt(tr[W*noiseCov*W]), e.g. for constructing
%                             pseudo-z-stats (see Vrba and Robinson 2001). 
%                             'none' applies no normalisation
%        .Normalise.leadFields - [true] or false normalises lead fields as
%                                a compensation for depth bias. 
%
%        .gammaPrior - form of prior on source variance parameters. Can
%                      be 'uniform-sd   - flat on sqrt(gamma) [Default].
%                         'uniform-variance' - flat on gamma.
%                         'log-uniform' - flat on log-gamma (this is the 
%                                         Jeffrey's prior but produces an 
%                                         improper posterior).
%                         'cauchy'      - weakly informative prior with a 
%                                         peak at zero, flat centre and
%                                         tail-off.
%                         'normal'      - on gamma
%                         'log-normal'  - on gamma
%                      or a handle to a function 
%                        lp(gamma) = sum_i (log p(gamma_i))
%                      which returns the log of a normalised prior density 
%                      for any scalar or vector gamma >= 0.
%


% References:
%  Wipf and Nagarajan. A unified Bayesian framework for MEG/EEG source imaging. Neuroimage (2009) vol. 44 (3) pp. 947-66
%
%   Dale, A.M. & Sereno, M.I.(1993) "Improved localization of cortical
%   activity by combining EEG and MEG with MRI cortical surface
%   reconstruction: a linear approach," J. Cogn. Neurosci 5, pp. 162--176.
%
%   Pascual-Marqui, R.D. (2002) "Standardized low-resolution brain
%   electromagnetic tomography (sLORETA): technical details," Methods Find.
%   Exp. Clin. Pharmacol. 24 (Suppl D) pp. 5--12.
%
%   Hamalainen, M.S., Lin, F. & Mosher, J.C. (2010) "Anatomically and
%   functionally constrained minimum-norm estimates." In Hansen, P.C.,
%   Kringelbach, M.L. & Salmelin, R. (ed.) "MEG: an introduction to
%   methods," Oxford University Press, Oxford, UK, pp. 186--215.
%
%   Wipf and Nagarajan (2007). Beamforming using the relevance vector
%   machine. Proc. 24th Int. Conf. on Machine Learning
%
%   Vrba and Robinson (2001). Signal processing in magnetoencephalography.
%   Methods 25(2), 249--271. 
%
%   Wens, V. et al. (2015) "A geometric correction scheme for spatial leakage effects in MEG/EEG seed-based functional connectivity mapping," Hum. Brain. Mapp. (In review)


%   Copyright 2015 OHBA
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


%   $LastChangedBy$
%   $Revision$
%   $LastChangedDate$
%   Contact: giles.colclough@magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 19-Jan-2015 14:40:17

global DEBUG
DEBUG = true;

%% Input Checking
% Sensor data covariance 
SensorData.nSensors = ROInets.rows(SensorData.cov);
validateattributes(SensorData.cov, {'numeric'},                 ...
                   {'2d', 'real', 'ncols', SensorData.nSensors}, ...
                   mfilename, 'SensorData.cov', 1);
% Check for positive definiteness? No - sensor data may not be full rank
% assert(ROInets.isposdef(SensorData.cov),              ...
%        [mfilename ':SensorDataCovarianceNotPosDef'], ...
%        'Data covariance matrix SensorData.cov should be positive definite. \n');

% lead fields
D = 3; % dimensions in space
% at the moment, expecting 3D lead fields only
assert(isequal(LeadFields.nDims, D),                 ...
       [mfilename ':LeadFieldsWrongDimensionality'], ...
       'Lead field dimensionality should be for 3D space. \n');
   
% check matrix size matches dimension fields and correct number of sensors
validateattributes(LeadFields.lf, {'numeric'},                        ...
                   {'2d', 'real', 'nrows', SensorData.nSensors,       ...
                    'ncols', LeadFields.nSources * LeadFields.nDims}, ...
                   mfilename, 'LeadFields.lf');

% set default options
Options = parse_options(Options);

%% Parse noise matrix
% noise covariance can be specified as the identity, the diagonal of the 
% data covariance, or passed in as a measured or estimated quantity. 
% For the first two cases, a global scaling is possible with parameter
% lambda. 
Noise = parse_noise(Noise, diag(SensorData.cov));

%% Lead field normalisation
% normalisation of the leadfields is often performed to compensate to some
% extent for the depth bias. 
LeadFields = normalise_leadfields(LeadFields, Options);

    
%% Source model and weights
fprintf('%s: computing weights using model %s. \n', ...
        mfilename, Options.sourceModel);
    
switch lower(Options.sourceModel)
    case 'wens'
        W_3d = wens_estimate(SensorData, Noise, LeadFields);
    case 'mne'
        W_3d = mne_estimate(SensorData, Noise, LeadFields, ...
                            Options.gammaPrior.fn);
    case 'mne-ina'
        W_3d = mne_estimate_ina(SensorData, Noise, LeadFields, ...
                                Options.gammaPrior.fn);
    case 'mne-scaled-noise'
        W_3d = mne_estimate_scale_noise(SensorData, Noise, LeadFields, ...
                                        Options.gammaPrior.fn);
    case 'sparse-bayes'
        W_3d = sparse_bayes_estimate(SensorData, Noise, LeadFields, ...
                                     Options.gammaPrior);
    case 'rvm-beamformer'
        W_3d = rvm_beamformer(SensorData, Noise, LeadFields, ...
                              Options.gammaPrior);
    otherwise
        error([mfilename ':InvalidSourceMethod'], ...
              'The chosen source method %s is not recognised. \n', ...
              Options.method);
end%switch

   
%% Parse into cell array and normalise weights
for iVox = LeadFields.nSources:-1:1, % initialise matrices by looping backwards
    
    % Partition weights for this dipole
    dipInd = (D*iVox - (D-1)):(D*iVox); % relevant indices for this dipole. LeadFields.nDims and D should be identical (checked above)
    Ws     = W_3d(dipInd,:);            % single-dipole weights
    
    if Options.ReduceRank.weights,
        % project onto direction of maximum variance using PCA. 
        dipCov         = Ws * SensorData.cov * Ws.';   % single-dipole covariance 3x3
        [ns,~]         = eigs(dipCov, [], 1);          % ns is direction of maximum source variance. This code seems about as fast as [ns,~] = fast_svds(dipCov, 1);
        W_nonorm{iVox} = ns.' * Ws;                    % weights     for this dipole projected onto direction of maximum variance
        lf{iVox}       = LeadFields.lf(:,dipInd) * ns; % lead fields for this dipole projected onto direction of maximum variance
    else
        W_nonorm{iVox} = Ws;
        lf{iVox}       = LeadFields.lf(:,dipInd);
    end%if
    
    % apply weights normalisation
    switch lower(Options.Normalise.weights)
        case 'sloreta'
            % normalising constant for depth bias using sLORETA
            % Wens et al. sec 2.4
            lambda_s = sqrt(trace(W_nonorm{iVox}                       ...
                                  * empirical_bayes_cov(Noise.cov,     ...
                                                        sourceCov,     ...
                                                        LeadFields.lf) ...
                                  * W_nonorm{iVox}.')); 
        case 'norm'
            % Normalise by norm of weights
            lambda_s = norm(W_nonorm{iVox});
            
        case 'noise-proj'
            % Normalise by projected noise power. 
            lambda_s = sqrt(trace(W_nonorm{iVox} * Noise.cov * W_nonorm{iVox}.'));
        case 'none'
            % apply no normalisation
            lambda_s = 1.0;
        otherwise
            error([mfilename ':InvalidNormalisationMethod'], ...
                  'Normalisation method %s is invalid. \n',  ...
                  Options.Normalise);
    end%switch
    
    W{iVox} = W_nonorm{iVox} ./ lambda_s; % normalised scalar weights vector for this dipole
end%loop over voxels

end%osl_inverse_mne_weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sensorCov = empirical_bayes_cov(noiseCov, sourceCov, lf3d)
%EMPIRICAL_BAYES_COV regularised estimate of data covariance in sensors
%
% SENSORCOV = EMPIRICAL_BAYES_COV(NOISECOV, SOURCECOV, 3D_LEADFIELDS)
%   returns the empirical bayes estimate of the sensor covariances based on
%   the noise covariance matrix, the modelled source covariance matrix and
%   the lead field matrix in 3D. 

sensorCov = noiseCov + lf3d * sourceCov * lf3d.';

% ensure real - sometimes goes a bit off. 
% sensorCov = real(sensorCov);
sensorCov = (sensorCov + sensorCov')./2.0;

end%empirical_bayes_cov





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W_3d = estimate_weights(noiseCov, sourceCov, lf3d, rho)
%ESTIMATE_WEIGHTS for reconstruction of 3d sources
%
% W = ESTIMATE_WEIGHTS(NOISECOV, SOURCECOV, 3D_LEADFIELDS) estimates
%   weights W for estimating sources S = W*B in 3d. 
%
% W = ESTIMATE_WEIGHTS(NOISECOV, SOURCECOV, 3D_LEADFIELDS, RHO) estimates
%   weights W for estimating sources S = W*B in 3d when there is an
%   additional scaling factor rho between the measured noise and measured
%   data.
%
% equation 13 from Wipf and Nagarajan (2009).
% given B = (LS + e)./rho, and S ~ N(0,sourceCov) then S is normally
% distributed with mean 
%    W * B = rho sourceCov L' (noiseCov + L sourceCov L')^{-1} B

% set default rho as identity
if nargin < 4 || ~exist('rho', 'var') || isempty(rho),
    rho = 1;
end%if

W_3d = rho .* sourceCov * lf3d.' * ...      
       osl_cholinv(empirical_bayes_cov(noiseCov, sourceCov, lf3d));
   
end%estimate_weights





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W_3d = beamformer_weights(noiseCov, sourceCov, lf3d, rho)
%ESTIMATE_WEIGHTS for reconstruction of 3d sources
%
% W = BEAMFORMER_WEIGHTS(NOISECOV, SOURCECOV, 3D_LEADFIELDS) estimates
%   weights W for estimating sources S = W*B in 3d. 
%
% W = BEAMFORMER_WEIGHTS(NOISECOV, SOURCECOV, 3D_LEADFIELDS, RHO) estimates
%   weights W for estimating sources S = W*B in 3d when there is an
%   additional scaling factor rho between the measured noise and measured
%   data.
%
% Vrba and Robinson 2001

% set default rho as identity
if nargin < 4 || ~exist('rho', 'var') || isempty(rho),
    rho = 1;
end%if

invDataCov = osl_cholinv(empirical_bayes_cov(noiseCov, sourceCov, lf3d));

W_3d = rho .* lf3d.' * invDataCov ./ trace(lf3d.' * invDataCov * lf3d);   
end%estimate_weights





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = wens_estimate(SensorData, Noise, LeadFields)
%WENS_ESTIMATE reconstruction weights by Vincent Wens' method
% regularization parameter from Wens et al. Sec 2.4

% use property sum(eig(B, A)) = trace(inv(A) * B)
k = sum(eig(LeadFields.lf * LeadFields.lf.', Noise.cov)) ./ ...
     (sum(eig(SensorData.cov, Noise.cov)) - SensorData.nSensors); 

% In our formulation, gamma = 1/k and C = I
sourceCov = (1.0 ./ real(k)) * speye(LeadFields.nDims * LeadFields.nSources);

% Extract weights for 3d sources
W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf);

end%wens_estimate





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W, logGamma] = mne_estimate(Data, Noise, LeadFields, prior)
%MNE_ESIMATE regularized source reconstruction weights 
%
% Estimates source covariance for MNE, using a white covariance matrix,
% with a single scaling parameter. 
% Use Eq 18 from Wipf and Nagarajan (2009), optimising w.r.t scalar gamma
global DEBUG

sourceCovFn     = @(logGamma) exp(real(logGamma)) .* ...
                              speye(LeadFields.nDims * LeadFields.nSources);
                                    
% variance of data can change by factor of e^20 relative to noise in search
% space
% note that for priors flat in log-space, which create improper posteriors,
% the inference can be quite sensitive to the lower bound used here.
noiseToLFScale = mean(log(diag(Noise.cov))) - ...
                 mean(log(sum(LeadFields.lf.^2)));
logGammaBound  = noiseToLFScale + [-20 50]; 

% if the prior has a scale, set it to be about e^7 above noise
scale = exp(noiseToLFScale + 7);

optimise_target = @(logGamma) optimise_target_single_prior(Data, Noise, ...
                          LeadFields, logGamma, prior, scale, sourceCovFn);
                      
% optimise using a golden section search and parabolic interpolation
logGamma      = fminbnd(optimise_target, logGammaBound(1), logGammaBound(2));
sourceCov     = sourceCovFn(logGamma);

% Extract weights for 3d sources
W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf);


if DEBUG,
    % plot
    lg = log(logspace(log10(exp(logGammaBound(1))), log10(exp(logGammaBound(2))), 50));
    L  = arrayfun(optimise_target, lg);
    figure('Color', 'w', 'Name', 'Optimisation target for log(gamma)');
    semilogy(lg, L);
    xlabel('Log (\gamma)');
    ylabel('L');    
    
    % display result
    fprintf('\nBounds on logGamma: %0.2g, %0.2g. \n', logGammaBound);
    fprintf('Chosen logGamma: %0.6g. \n', logGamma);
    
end%if DEBUG
end%mne_Jeffreys_prior_estimate








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = mne_estimate_scale_noise(Data, Noise, LeadFields, priorFn)
%MNE_DOUBLE_SEARCH_ESIMATE regularized source estimates with noise scaling
%
% Estimates source covariance for MNE, using a white covariance matrix,
% with a single scaling parameter. 
% Allows an additional parameter to control scaling between noise and data
% as rho*B = LS + e
global DEBUG

% It seems as though this functionality will not work well. No minimum for
% rho, unless forced by prior?
warning([mfilename ':scale_noise'],                    ...
        ['The scale-noise mne estimate is unstable. ', ...
         'It may not find a minimum for rho. \n']);

sourceCovFn     = @(logGamma) exp(real(logGamma)) .* ...
                              speye(LeadFields.nDims * LeadFields.nSources);
                                    
% put noise and sources are on same scale
noiseToLFScale = mean(log(diag(Noise.cov))) - ...
                 mean(log(sum(LeadFields.lf.^2)));

% put data and noise on scale set by mean eigenvalue
dataEigVals      = eig(Data.cov);
noiseEigVals     = eig(Noise.cov);
NoiseToDataScale = log(mean(noiseEigVals)) - log(mean(dataEigVals));

% set initial values
logGammaInit = noiseToLFScale; 
logRhoInit   = NoiseToDataScale./2;
paramsInit   = [logGammaInit; logRhoInit];

% priors forms are the same for rho and gamma
priors = {priorFn; priorFn};

% If there is a scale on the priors, we want to constrain the space to
% moderately sensible values
% If we reckon our initial guesses are any good, let's constrain the
% variance parameters to be within a factor of exp(7).
scales = exp(paramsInit + 7);

optimise_target = @(params) optimise_target_double_prior(Data, Noise,        ...
                                                         LeadFields, params, ...
                                                         priors, scales,     ...
                                                         sourceCovFn);

% set optimisation params
maxIter = 4000;
                                                     
% optimise using a multivariate nonlinear Nelder-Mead minimization
[params, fval, success] = fminsearch(optimise_target, paramsInit,     ...
                                     optimset('MaxFunEvals', maxIter, ...
                                              'MaxIter', maxIter));

logGamma  = params(1);
sourceCov = sourceCovFn(logGamma);
rho       = exp(params(2));

% Extract weights for 3d sources
W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf, rho);

if DEBUG || ~success,
    % use a grid search
    lrBound = real(log10(rho)) + [-7 7];
    lgBound = real(log10(exp(logGamma))) + [-10 10];
    nR      = 20;
    nG      = 30;
    lr      = log(logspace(lrBound(1), lrBound(2), nR));
    lg      = log(logspace(lgBound(1), lgBound(2), nG));
    for iR = length(lr):-1:1,
        for iG = length(lg):-1:1,
            L(iG,iR) = optimise_target([lg(iG), lr(iR)]);
        end%for
    end%for
    
    figure('Color', 'w', 'Name', 'Optimisation target for log(gamma)');
    surfl(real(lr), real(lg), real(L));
    colormap pink
    xlabel('Log (\rho)');
    ylabel('Log (\gamma)');
    zlabel('L');
    
    % best gamma and rho?
    [searchFval, iMin] = min(L(:));
    if searchFval < fval,
        fprintf('%s: Using best parameters from grid search. \n', mfilename);
        [iG,iR]   = ind2sub(size(L), iMin);
        logGamma  = lg(iG);
        sourceCov = sourceCovFn(logGamma);
        rho       = exp(lr(iR));
        
        W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf, rho);
    end%if
end%if DEBUG
end%mne_estimate_scale_noise







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = mne_estimate_ina(Data, Noise, LeadFields, prior)
%MNE_ESTIMATE_INA regularized source estimates
%
% Estimates source distribution using mean of marginal distribution on
% sources. Marginal distribution is approximated by averaging over several
% discrete points in gamma-space, rather than taking the MAP estimate of
% gamma. This is the integrated nested approximation. 
% Returns weight W to estimate sources as S = W*B. 

global DEBUG

% Find or estimate peak of marginal gamma distribution
% either:
[~, logGammaPeak] = mne_estimate(Data, Noise, LeadFields, prior);
gammaPeak         = exp(logGammaPeak);
% or:
noiseToLFScale = mean(log(diag(Noise.cov))) - ...
                 mean(log(sum(LeadFields.lf.^2)));
% gammaPeak      = exp(noiseToLFScale);

% define inverse weights function 
nSourceElements = LeadFields.nDims * LeadFields.nSources;
weights         = @(g) estimate_weights(Noise.cov,                 ...
                                        g.*speye(nSourceElements), ...
                                        LeadFields.lf);

% if the prior has a scale, set it to be about e^7 above noise
scale       = exp(noiseToLFScale + 7);
sourceCovFn = @(logGamma) exp(real(logGamma)) .* speye(nSourceElements);

% integrate with clenshaw-curtis integration in 3 bands
nLower = 2^7  + 1;
nInner = 2^12 + 1;
nUpper = 2^6  + 1;
gammaBounds = [gammaPeak./1e40, gammaPeak./1e9, gammaPeak*1e6, gammaPeak*1e30];
[g1, w1] = clen_curt_points(nLower, gammaBounds(1), gammaBounds(2)*0.9999);
[g2, w2] = clen_curt_points(nInner, gammaBounds(2), gammaBounds(3)*0.9999);
[g3, w3] = clen_curt_points(nUpper, gammaBounds(3), gammaBounds(4));
g      = cat(1, g1, g2, g3);
w      = cat(1, w1, w2, w3);
[g,is] = sort(g);
w      = w(is);

% we want to avoid exact zeros in x when working in log-space
exactZeros = 0 == g;
g(exactZeros) = [];
w(exactZeros) = [];
N = length(g);

% set up tracking plot
if DEBUG,
    fprintf('Finding values of logp(g|B). \n');
    flowBound = eps(min(g));
    logPrior = arrayfun(@(lx) prior(lx, scale), log(g+flowBound));
    
    figure('Name', 'p(g)', 'Color', 'w');
    plot(log(g), exp(logPrior - log_sum_exp(logPrior)), 'r');
    xlabel('Log(\gamma)');
    hold on;
end%if

% declare memory
W      = zeros(nSourceElements, Data.nSensors);
log_pg = zeros(N,1);

% find coefficients for expected weights at each x
for i = 1:N,
    if ~mod(i, 100), 
        fprintf('%s: marginal contribution %d of %d. \n', ...
                mfilename, i, N);
    end%if
    log_pg(i) = -0.5*optimise_target_single_prior(Data, Noise, ...
                          LeadFields, log(g(i)), prior, scale, sourceCovFn); %logp_g_on_B(x(i));
end%for
log_wpg    = log_pg + log(w);
log_Exp_pg = log_sum_exp(log_wpg);

if DEBUG,
    figure('Name', 'p(g|B)', 'Color', 'w');
    plot(log(g), log_pg - log_Exp_pg, 'k');
    xlabel('Log(\gamma)');
end%if

% loop over integration points, accumulating result
for iInt = 1:N,
    W = W + exp(log_wpg(iInt) - log_Exp_pg) .* weights(g(iInt));
    
    % write out progress
    if DEBUG 
        if ~mod(iInt, 100), 
            fprintf('%s: INA contribution %d of %d. \n', ...
                    mfilename, iInt, N);
        end%if
    end%if
end%for

if DEBUG,
    fprintf('%s: normalisation check: sum of weighting factors: %0.3g. \n', ...
            mfilename, exp(log_sum_exp(log_wpg - log_Exp_pg)));
        
    gammaMapW             = weights(gammaPeak);
    fractionalImprovement = norm(W - gammaMapW)./norm(W) * 100;
    
    gammaMean = sum(exp(log_wpg - log_Exp_pg) .* g);
        
    fprintf('%s: difference between INA estimate and gamma-MAP: %0.2g%%. \n', ...
            mfilename, fractionalImprovement);
    fprintf(['%s: gamma-MAP = %0.4g,\tE[gamma] = %0.4g,', ...
             '\tdifference: %0.3g%%. \n'],                ...
            mfilename, gammaPeak, gammaMean,              ...
            (gammaMean - gammaPeak)./gammaMean * 100);
end%if
end%mne_estimate_ina
%{
% % Previously used a nested function:
%     function logp = logp_g_on_B(g)
%         %LOGP_G_ON_B
%         weightsCov = weights(g)'*weights(g);
%         IminusLW = (speye(Data.nSensors) - LeadFields.lf*weights(g));
%         logp     = - 0.5 * Data.nSamples * Data.nSensors * log(2*pi)                          ...
%                    - 0.5 * Data.nSamples * logDetNoise                                        ...
%                    - 0.5 * trace(IminusLW * Data.cov * IminusLW.' * noiseInv) * Data.nSamples ...
%                    - 0.5 * Data.nSamples * nSourceElements * log(2*pi*g)                      ...
%                    - 0.5 * sum(sum(weightsCov .* Data.cov)) * Data.nSamples ./ g              ...
%                    + prior(log(g),scale);
%     end%logp_g_on_B
%}






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = sparse_bayes_estimate(Data, Noise, LeadFields, GammaPrior)
%SPARSE_BAYES_ESIMATE weights for sparse bayes source reconstruction
%
% Estimates source covariance using sparse Bayes estimate: allow for
% source variance at each location and in each orientation in space. 
%
% To use this algorithm, -log p(gamma) MUST BE CONCAVE for each gamma.
%
% Ref: Wipf and Nagarajan, 2009

% estimate regularized covariance
sourceCov = sparse_bayes_covariance(Data, Noise, LeadFields, GammaPrior);

% Extract weights for 3d sources
W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf);
end%sparse_bayes_estimate






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = rvm_beamformer(Data, Noise, LeadFields, GammaPrior)
%RVM_BEAMFORMER weights for relevance vector machine beamforming
%
% Estimates source covariance using sparse Bayes estimate: allow for
% source variance at each location and in each orientation in space,
% according to Wipf and Nagarajan 2007, but updated to use a full noise
% covariance estimate. 
%
% MAP optimisation of source variances occurs using formulae from Wipf and
% Nagarajan, 2009. 
%
% To use this algorithm, -log p(gamma) MUST BE CONCAVE for each gamma.

% estimate regularized covariance
% Wipf and Nagarajan 2007 eq 10
sourceCov = sparse_bayes_covariance(Data, Noise, LeadFields, GammaPrior);

% Extract weights for 3d sources
% Wipf and Nagarajan 2007 Eq 7
W = beamformer_weights(Noise.cov, sourceCov, LeadFields.lf);
end%rvm_beamformer








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sourceCov = sparse_bayes_covariance(Data, Noise, LeadFields, ...
                                             GammaPrior)
%SPARSE_BAYES_COVARIANCE regularized source covariance 
%
% Estimates source covariance using sparse Bayes estimate: allow for
% source variance at each location and in each orientation in space. 
%
% To use this algorithm, -log p(gamma) MUST BE CONCAVE for each gamma.
%
% Ref: Wipf and Nagarajan, 2009
global DEBUG
nSourceElements = LeadFields.nDims * LeadFields.nSources;
sourceCovFn     = @(gamma) spdiags(real(gamma), ...
                                   0, nSourceElements, nSourceElements);
logGammaInit    = (mean(log(diag(Noise.cov)))          ...
                   - mean(log(sum(LeadFields.lf.^2)))) ...
                  .* ones(nSourceElements, 1);
              
% set scale for priors on gamma
scale = exp(logGammaInit + 5);
              
% use updates from gamma-MAP in Wipf and Nagarajan (2009)
oldGamma   = exp(logGammaInit);
nIterMax   = 1000;
iIter      = 0;
deltaGamma = NaN(nIterMax,1);
gamma      = zeros(size(oldGamma));
L          = NaN(nIterMax,1);

HALT_CONDITION = 1e-4; %fractional change in norm(gamma) to class as convergence

while iIter <= nIterMax,
    % some monitoring behaviour
    if DEBUG,
        % monitor cost function Eq 18 in Wipf and Nagarajan
        if iIter > 0,
            L(iIter) = optimise_target_single_prior(Data, Noise, LeadFields,           ...
                                                    log(oldGamma), GammaPrior.fn,      ...
                                                    scale, @(lg) sourceCovFn(exp(lg)), ...
                                                    true);
        end%if
                                   
        if (~mod(iIter, 5) || 1 == iIter) && 0 ~= iIter,
        fprintf('%s: SB g-MAP iteration no %4d, deltaGamma %0.5g, L %0.6g.\n', mfilename, iIter, deltaGamma(iIter), L(iIter));
        end%if
        if  0 == iIter,
            ff   = figure('Color', 'w', 'Name', 'Updating gamma estimates');
            rmFF = onCleanup(@() close(ff));
            hist(log(gamma), 3000);
            xlabel('Log(\gamma)');
        elseif ~mod(iIter, 5) || (1 == iIter),
            nzg = gamma~=0;
            hist(log(gamma(nzg)), 3000);
            title(sprintf('fraction non-zero: %0.2f%%', ...
                          sum(nzg)./nSourceElements * 100));
            drawnow;
        end%if
    end%if
        
    % start of loop proper
    iIter = iIter + 1;
    
    % use Factorize object to hold inverse of regularized covariance
    % without computation
    invSigmaB = osl_cholinv(empirical_bayes_cov(Noise.cov, sourceCovFn(oldGamma), LeadFields.lf));
    
    % update rules
    % Wipf and Nagarajan 2009 Eq. 31
    for i = nSourceElements:-1:1,
        gamma(i) = oldGamma(i)                                                                                   ...
                   *  sqrt(trace(LeadFields.lf(:,i).' * invSigmaB * Data.cov * invSigmaB *  LeadFields.lf(:,i))) ... % Frobenius norm term
                   ./ sqrt(trace(LeadFields.lf(:,i).' * invSigmaB * LeadFields.lf(:,i)                           ...
                                 - GammaPrior.diff_fn(oldGamma(i), scale) ./ Data.nSamples));                        % Final trace term
    end%for 
    
    % monitor change
    deltaGamma(iIter) = norm(gamma - oldGamma) ./ norm(oldGamma);
    
    % test for convergence
    isConverged = deltaGamma(iIter) < HALT_CONDITION;
    if isConverged,
        break
    end%if
    
    % re-assign old gamma
    oldGamma = gamma;
end%while

% tidy up
deltaGamma(isnan(deltaGamma)) = [];

% check if hit max iter
if iIter >= nIterMax,
    warning([mfilename ':gMAPMaxIter'], ...
            'Max number of iterations, %d, reached without convergence. \n', ...
            nIterMax);
end%if

% final result
sourceCov = sourceCovFn(gamma);


if DEBUG
    delete(rmFF); % delete old figure
    
    figure('Name', 'Sparse-Bayes convergence', 'Color', 'w');
    semilogy(deltaGamma, 'r', 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', 14);
    ylabel('Fractional change in norm(\gamma)', 'FontSize', 14);
    
    figure('Name', 'Gamma estimates', 'Color', 'w');
    width  = 0.8;
    nzg = gamma~=0;
    [n,x] = hist(log(gamma(nzg)), 3000);
    hh    = bar(x, n, width);
    set(get(hh, 'Children'),              ...
        'FaceColor', [120, 120, 120]/255, ... % set to be grey
        'FaceAlpha', 0.4,                 ...
        'EdgeColor', 'none');
    title(sprintf('fraction non-zero: %0.2f%%', ...
                  sum(nzg)./nSourceElements * 100));
    xlabel('Log(\gamma)');
end%if
end%sparse_bayes_covariance





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = optimise_target_single_prior(Data, Noise, LeadFields, ...
                                          logGamma, prior, scale,  ...
                                          sourceCovFn, quiet)
%OPTIMISE_TARGET_single_PRIOR
% minimise L to find gamma-MAP solution
global DEBUG
persistent iTARGETCALL
if nargin < 8 || ~exist('quiet', 'var'),
    quiet = false;
end%if

Sigma_EB = empirical_bayes_cov(Noise.cov, sourceCovFn(logGamma), ...
                               LeadFields.lf);
                           
% we want to maximise
% log p(g|B) = -0.5 Tr(BB' Sigma_b^{-1}) - n/2 logdet Sigma_b + log p(g).

try
    logDetSigma = ROInets.logdet(Sigma_EB, 'chol');
    
catch ME
    % some errors with non pos def matrices occuring, which is surprising. 
    if strcmp(ME.identifier, 'MATLAB:posdef'),
        warning([mfilename ':OptimTarget:posdef'], ...
                'Regularised covariance not positive defninite. \n');
        
        logDetSigma = ROInets.logdet(Sigma_EB);        
    else
        rethrow(ME);
    end%if
end%try

% use property sum(eig(B, A)) = trace(inv(A) * B)
% or trace(AB) = sum(sum(A .* B')) (and covariance matrices are symmetric)
L = real((trace(Data.cov * osl_cholinv(Sigma_EB)) ...                  % faster than elementwise product or sum(eig()). 
     + logDetSigma) * Data.nSamples - 2 * prior(logGamma, scale));       


if DEBUG && ~quiet,
    if ~exist('iTARGETCALL', 'var') || isempty(iTARGETCALL),
        iTARGETCALL = 1;
    else
        iTARGETCALL = iTARGETCALL + 1;
    end%if
    fprintf(['Call to optim fn %4.0d: L = %0.8G, logGamma = %0.4G, ', ...
             'logdet(Sigma_EB) = %0.6G. \n'],                         ...
            iTARGETCALL, L, logGamma, logDetSigma);
end%if DEBUG
end%optimise_target_single_prior





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = optimise_target_double_prior(Data, Noise, LeadFields, params, ...
                                          priors, priorScales, sourceCovFn)
%OPTIMISE_TARGET_DOUBLE_PRIOR
% minimise L to find gamma-MAP solution for S under B = (LS + e)/r.

global     DEBUG
persistent iTARGETCALL

logGamma = real(params(1));
logRho   = real(params(2));
pg       = priors{1};
pr       = priors{2};

Sigma_EB = empirical_bayes_cov(Noise.cov, sourceCovFn(logGamma), ...
                               LeadFields.lf);

try
    logDetSigma = ROInets.logdet(Sigma_EB, 'chol');
    
catch ME
    % some errors with non pos def matrices occuring, which is surprising. 
    if strcmp(ME.identifier, 'MATLAB:posdef'),
        warning([mfilename ':OptimTarget:posdef'], ...
                'Regularised covariance not positive defninite. \n');
        
        logDetSigma = ROInets.logdet(Sigma_EB);        
    else
        rethrow(ME);
    end%if
end%try
                           
% we want to maximise
% log p(r,g|B) = -0.5 r^2 Tr(BB' Sigma_b^{-1}) - n/2 logdet Sigma_b  
%                + n*log(r) + log p(r) + log p(g).
%
% use property sum(eig(B, A)) = trace(inv(A) * B)
% or trace(AB) = sum(sum(A .* B')) (and covariance matrices are symmetric)
L = (exp(2*logRho) * trace(Data.cov * osl_cholinv(Sigma_EB)) ...  % faster than elementwise product or sum(eig()). 
     + logDetSigma - 2*logRho) .* Data.nSamples  ...
    - 2*pg(logGamma, priorScales(1)) - 2*pr(logRho, priorScales(2)); 


if DEBUG,
    if ~exist('iTARGETCALL', 'var') || isempty(iTARGETCALL),
        iTARGETCALL = 1;
    else
        iTARGETCALL = iTARGETCALL + 1;
    end%if
    fprintf(['Call to optim fn %4.0d: L = %0.8G, logGamma = %0.4G, ', ...
             'logRho = %0.4G, logdet(Sigma_EB) = %0.6G. \n'],         ...
            iTARGETCALL, L, logGamma, logRho,                         ...
            logDetSigma);
end%if DEBUG
end%optimise_target_double_prior





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lp, diff_lp] = parse_priors(priorChoice)
%PARSE_PRIORS creates function handle for priors on hierarchical variance
%parameters. 
%
% Some priors require scales, so lp takes two arguments, lp(loggamma, scale).
% May also require d(lp)/dg. (= - df(g)/dg in Wipf and Nagarajan)

% pass on any passed-in functions
if isa(priorChoice, 'function_handle'),
    lp      = @(log_g,~) priorChoice(exp(log_g));
    diff_lp = @(~,~) error([mfilename ':NoDiffLP'], ...
                           'No differential of gamma prior passed in.\n');
    return
end%if

switch lower(priorChoice)
    case 'uniform-variance'
        % p(g) ~ 1
        lp      = @(log_g,~) 0; % ignore additive constant
        diff_lp = @(g,~)     0;
        
    case 'uniform-sd'
        %p(sqrt(g)) ~ 1 --> p(g) ~ 0.5 / sqrt(g)
        lp      = @(log_g,~) - 0.5 * sum(log_g); % ignore additive constant
        diff_lp = @(g,~)     - 0.5 ./ sum(g);
        
    case 'uniform-log'
        % p(g) ~ 1/g
        lp      = @(log_g,~) - sum(log_g); % ignore additive constants as improper unless bounded
        diff_lp = @(g,~)     - 1.0 ./ sum(g);
        
    case 'cauchy'
        % p(sqrt(g) ~ 1/pi sqrt(A) / (g+A) 
        % --> p(g) ~ 1/(2 pi) sqrt(A) / ((A+g) sqrt(g))
        lp      = @(log_g,A) sum(-log(2*pi) + 0.5*log(A) ...
                                 - 0.5*log_g - log(exp(log_g) + A));
                             
        diff_lp = @(g,A)     sum(-0.5 ./ g - 1.0 ./ (g + A));
        
    case 'log-normal' % scale on same space as g, so s.d. is log(A).
        lp      = @(log_g,A) - 0.5 * (log_g.' * log_g) ./ log(A)^2 ...
                             - length(log_g) * 0.5 * log(2*pi*log(A)^2); 
                         
        diff_lp = @(g,A)     sum( - log(g) ./ (log(A)^2 * g));   
        
    case 'normal'
        lp      = @(log_g,A) - 0.5 * (exp(log_g).' * exp(log_g)) ./ A^2 ...
                             - length(log_g) * 0.5 * log(2*pi*A^2);
                    
        diff_lp = @(g,A)     - sum(g) ./ A^2;
                    
    otherwise
        error([mfilename ':UnexpectedPrior'],                    ...
              'The choice of prior, %s, was not recognised. \n', ...
              priorChoice);
end%switch
end%parse_priors





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LeadFields = normalise_leadfields(LeadFields, Options)
%NORMALISE_LEADFIELDS
%
% LEADFIELDS = NORMALISE_LEADFIELDS(LEADFIELDS, OPTIONS) normalises lead
% fields if Options.Normalise.leadFields is true. 

if Options.Normalise.leadFields,
    LeadFields.isNormalised = true;
    
    % loop over sources and normalise each lead field
    for iSource = LeadFields.nSources:-1:1,
        % extract relevant lead field
        sourceXYZ = ((iSource-1) * LeadFields.nDims + 1):(iSource * LeadFields.nDims);
        lf        = LeadFields.lf(:,sourceXYZ);
        
        % take magnitude
        LeadFields.lfMag(iSource)  = norm(lf);
        
        % normalise
        LeadFields.lf(:,sourceXYZ) = lf ./ LeadFields.lfMag(iSource);
    end%for
else
    LeadFields.isNormalised = false;
end%if
end%normalise_leadfields





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Noise = parse_noise(Noise, diagDataCov)
%PARSE_NOISE sets noise covariance matrix
%
% NOISE = PARSE_NOISE(NOISE, DIAGONAL_OF_DATA_COVARIANCE)
%
% noise covariance can be specified as the identity, the diagonal of the 
% data covariance, or passed in as a measured or estimated quantity. 
% For the first two cases, a global scaling is possible with parameter
% lambda. 

methodMixWarning = ['You have passed in a noise covariance matrix but ', ...
                    'not specified the empirical noise method. \n',      ...
                    'The noise covariance matrix will be ignored. \n'];
hasNoiseCov      = isfield(Noise, 'cov') && ~isempty(Noise.cov);

nSensors       = length(diagDataCov);

switch lower(Noise.model)
    case 'white'
        Noise.cov = Noise.lambda .* eye(nSensors);
        if hasNoiseCov,
            warning([mfilename ':noiseMethodMix'], methodMixWarning);
        end%if

    case 'diag_datacov'
        Noise.cov = Noise.lambda .* diag(diagDataCov);
        if hasNoiseCov,
            warning([mfilename ':noiseMethodMix'], methodMixWarning);
        end%if

    case 'empirical'
        % use Noise.cov passed in
        if ~hasNoiseCov,
            error([mfilename ':NoNoiseCov'],                        ...
                  ['Empirical noise method selected but no noise ', ...
                   'covariance provided. \n']);
        end%if
        validateattributes(Noise.cov, {'numeric'},                ...
                           {'real', 'size', [nSensors nSensors]}, ...
                           mfilename, 'Noise.cov');
        assert(ROInets.isposdef(Noise.cov), ...
               [mfilename ':NoiseCovNotPosDef'], ...
               'Noise.cov matrix must be positive definite. \n');

    otherwise % catch incorrect usage
        error([mfilename ':InvalidNoiseMethod'],     ...
              'Noise method %s not recognised. \n', ...
              Noise.method);
end%switch
end%parse_noise





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ParsedOptions = parse_options(Options)
%PARSE_OPTIONS sets default options
%
% OPTIONS = PARSE_OPTIONS(OPTIONS) manages settings for
% Options.ReduceRank.weights    - Boolean [1] converts to scalar weights
%                   .leadFields - Integer [3] dimensionality of lead fields
%        .sourceModel     - Wens (based on paper referenced below), 
%                           MNE (which uses gamma-MAP estimation), 
%                           MNE-INA (integrated nested approximation for
%                           averaging over the prior on gamma)
%                           MNE-scaled-noise (also estimates a parameter
%                           rho for the scaling between sensor data and
%                           noise)
%                           sparseBayes (gamma-MAP estimation of a
%                           sparseBayes solution)
%        .Normalise       - sLoreta, [norm], none
%        .gammaPrior - form of prior on source variance parameters. Can
%                      be 'uniform'     - flat on gamma [Default].
%                         'log-uniform' - flat on log-gamma (this is the 
%                                         Jeffrey's prior but produces an 
%                                         improper posterior).
%                         'Cauchy'      - weakly informative prior with a 
%                                         peak at zero, flat centre and
%                                         tail-off.
%                      or a handle to a function lp(gamma) which returns a
%                      log of a normalised prior density for any gamma>=0.

% set up input parser
IP               = inputParser;
IP.CaseSensitive = false;
IP.FunctionName  = mfilename;
IP.StructExpand  = true;  % If true, can pass parameter-value pairs in a struct
IP.KeepUnmatched = false; % If true, accept unexpected inputs

% declare allowed options
test = @(x) isstruct(x)         &&                         ...
       isfield(x, 'weights')    && islogical(x.weights) && ...
       isfield(x, 'leadFields') && isscalar(x.leadFields);
IP.addParamValue('ReduceRank', [], test);
test = @(x) ischar(validatestring(x,                                   ...
                                  {'Wens', 'MNE', 'MNE-INA',           ...
                                   'MNE-scaled-noise', 'sparse-Bayes', ...
                                   'rvm-beamformer'},                  ...
                                  mfilename, ...
                                  'Options.sourceModel'));
IP.addParamValue('sourceModel', 'MNE', test);
test = @(x) isstruct(x) &&                                                                     ...
            isfield(x, 'weights')    && ischar(validatestring(x.weights,                       ...
                                                              {'sLoreta', 'norm',              ...
                                                               'noise-proj', 'none'},          ...
                                                              mfilename,                       ...
                                                              'Options.Normalise.weights')) && ...
            isfield(x, 'leadFields') && islogical(x.leadFields);
        
IP.addParamValue('Normalise', 'norm', test);
test = @(x) isa(x, 'function_handle') || ...
            ischar(validatestring(x, {'uniform-sd', 'uniform-variance', ...
                                      'uniform-log', 'cauchy',          ...
                                      'normal', 'log-normal'},          ...
                                  mfilename, 'Options.gammaPrior'));
IP.addParamValue('gammaPrior', 'uniform-sd', test);

% parse inputs
IP.parse(Options);

ParsedOptions = IP.Results;

% set prior functions
[ParsedOptions.gammaPrior.fn, ...
 ParsedOptions.gammaPrior.diff_fn] = parse_priors(ParsedOptions.gammaPrior);
end%parse_options





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = clen_curt_points(N1,a,b)
%CLEN_CURT_POINTS Clenshaw Curtis quadrature
%
% [X, W] = CLEN_CURT_POINTS(N, A, B) computes points and weights for
%    clenshaw-curtis quadrature on the interval [A,B]. 1D integrals of f(x)
%    can be approximated using Sum_(i=1)^N w_i * f(x_i). 
%
% The algorithm is most efficient for choices of N = 2^p + 1
% J�rg Waldvogel, "Fast construction of the Fej�r and Clenshaw-Curtis 
% quadrature rules", BIT Numerical Mathematics 43 (1), p. 001-018 (2004).
%
% Adapted from code written by: Greg von Winckel - 02/12/2005
% Contact: gregvw(at)chtm(dot)unm(dot)edu


N   = N1-1; 
bma = b-a;
c   = zeros(N1,2);

c(1:2:N1, 1) = (2 ./ [1, 1 - (2:2:N).^2]).'; c(2,2)=1;

f = real(ifft([c(1:N1,:);c(N:-1:2,:)]));

w = bma * ([f(1,1); 2*f(2:N,1); f(N1,1)])./2;
x = 0.5 * ((b + a) + N * bma * f(1:N1,2));
end%clen_curt_points





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = log_sum_exp(x, dim)
%LOG_SUM_EXP avoids numerical underflow
%
% S = LOG_SUM_EXP(A) returns S = log(sum(exp(A)))
% 
% Default behaviour is to sum down columns. Use S = LOG_SUM_EXP(A,2) to sum
%   over rows instead. 

if nargin == 1,
    if(ROInets.rows(x) > 1)
      dim = 1;
    elseif(ROInets.cols(x) > 1)
      dim = 2;
    elseif isscalar(x)
        S = x; 
        return
    else
        dim = find(size(x), 1, 'first');
    end%if
end%if dim not provided

% find max in each column and subtract
maxVals = max(x(isfinite(x)), [], dim);
a       = bsxfun(@minus, x, maxVals);
S       = maxVals + log(sum(exp(a), dim));

% if elements of maxVals are infinite, use them in place
dodgyInds    = ~isfinite(maxVals);
S(dodgyInds) = maxVals(dodgyInds);

end%log_sum_exp
% [EOF]