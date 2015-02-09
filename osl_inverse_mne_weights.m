function [W, W_nonorm,lf] = osl_inverse_mne_weights(SensorData, LeadFields, Noise, Options)
%OSL_INVERSE_MNE_WEIGHTS creates weights for source reconstruction using MNE
%
% W = OSL_INVERSE_MNE_WEIGHTS(SensorData, LeadFields, NOISE, RESTRICTRANK)
%   calculates cell array
%
% SensorData.cov
%           .nSamples
%
% LeadFields.nSources
%           .nDims
%           .lf
% Noise.model
%      .lambda
%      .cov
%
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
%        .normalise       - sLoreta, [norm], none
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
% References:
%  Wipf and Nagarajan. A unified Bayesian framework for MEG/EEG source imaging. Neuroimage (2009) vol. 44 (3) pp. 947-66
%
%	Dale, A.M. & Sereno, M.I.(1993) "Improved localization of cortical
%	activity by combining EEG and MEG with MRI cortical surface
%	reconstruction: a linear approach," J. Cogn. Neurosci 5, pp. 162--176.
%
%	Pascual-Marqui, R.D. (2002) "Standardized low-resolution brain
%	electromagnetic tomography (sLORETA): technical details," Methods Find.
%	Exp. Clin. Pharmacol. 24 (Suppl D) pp. 5--12.
%
%   Hamalainen, M.S., Lin, F. & Mosher, J.C. (2010) "Anatomically and
%   functionally constrained minimum-norm estimates." In Hansen, P.C.,
%   Kringelbach, M.L. & Salmelin, R. (ed.) "MEG: an introduction to
%   methods," Oxford University Press, Oxford, UK, pp. 186--215.
%
%   Wens, V. et al. (2015) "A geometric correction scheme for spatial leakage effects in MEG/EEG seed-based functional connectivity mapping," Hum. Brain. Mapp. (In review)


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


%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 19-Jan-2015 14:40:17

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
    
%% Source model and weights
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
    case 'sparsebayes'
        W_3d = sparse_bayes_estimate(SensorData, Noise, LeadFields, ...
                                     Options.gammaPrior.fn);
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
    switch lower(Options.normalise)
        case 'sloreta'
            % normalising constant for depth bias using sLORETA
            % Wens et al. sec 2.4
            lambda_s = sqrt(W_nonorm{iVox}                              ...
                            * empirical_bayes_cov(Noise.cov, sourceCov, ...
                                                  LeadFields.lf)        ...
                            * W_nonorm{iVox}.'); 
        case 'norm'
            % normalise by norm of weights
            lambda_s = norm(W_nonorm{iVox});
        case 'none'
            % apply no normalisation
            lambda_s = 1.0;
        otherwise
            error([mfilename ':InvalidNormalisationMethod'], ...
                  'Normalisation method %s is invalid. \n',  ...
                  Options.normalise);
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
sensorCov = real(sensorCov);

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
       inverse(empirical_bayes_cov(noiseCov, sourceCov, lf3d));
   
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
% noiseToLFScale = median(log(diag(Noise.cov))) - ...
%                  median(log(sum(LeadFields.lf.^2)));
% gammaPeak      = exp(noiseToLFScale);

% define inverse weights function 
nSourceElements = LeadFields.nDims * LeadFields.nSources;
weights         = @(g) estimate_weights(Noise.cov,                 ...
                                        g.*speye(nSourceElements), ...
                                        LeadFields.lf);

% define p(gamma|B)
logDetNoise = logdet(Noise.cov, 'chol');
noiseInv    = inverse(Noise.cov);

p_g_on_B = @(g) exp(logp_g_on_B(g));

% integrate with clenshaw-curtis integration in 3 bands
nOuter = 2^14+1;
nInner = 2^19+1;
N      = nOuter*2 + nInner;
gammaBounds = [eps(gammaPeak), gammaPeak./1e5, gammaPeak*1e5, gammaPeak*1e50];
[x1, w1] = clen_curt_points(nOuter, gammaBounds(1), gammaBounds(2)*0.9999);
[x2, w2] = clen_curt_points(nInner, gammaBounds(2), gammaBounds(3)*0.9999);
[x3, w3] = clen_curt_points(nOuter, gammaBounds(3), gammaBounds(4));
x = cat(1, x1, x2, x3);
w = cat(1, w1, w2, w3);

Exp_W  = [];
Exp_pg = [];
for iInt = 1:N,
    if DEBUG && ~mod(iInt, 1000), 
        fprintf('%s: INA contribution %d of %d. \n', mfilename. iInt, N);
    end%if
    wpg = w(iInt) * p_g_on_B(x(iInt));
    Exp_W  = Exp_W  + wpg * weights(x(iInt));
    Exp_pg = Exp_pg + wpg;
end%for

W = Exp_W ./ Exp_pg;

    function logp = logp_g_on_B(g)
        %LOGP_G_ON_B
        weightsCov = weights(g)'*weights(g);
        IminusLW = (speye(nSourceElements) - LeadFields.lf*weights(g));
        logp     = - 0.5 * Data.nSamples * Data.nSensors * log(2*pi)                          ...
                   - 0.5 * Data.nSamples * logDetNoise                                        ...
                   - 0.5 * trace(IminusLW * Data.cov * IminusLW.' * noiseInv) * Data.nSamples ...
                   - 0.5 * Data.nSamples * nSourceElements * log(2*pi*g)                      ...
                   - 0.5 * sum(sum(weightsCov .* Data.cov)) * Data.nSamples ./ g              ...
                   + prior(log(g),scale);
    end%logp_g_on_B

end%mne_estimate_ina









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
noiseToLFScale = median(log(diag(Noise.cov))) - ...
                 median(log(sum(LeadFields.lf.^2)));
logGammaBound  = noiseToLFScale + [-10 20]; 

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
    lg = log(logspace(log10(exp(logGammaBound(1))), log10(exp(logGammaBound(2))), 50));
    L  = arrayfun(optimise_target, lg);
    figure('Color', 'w', 'Name', 'Optimisation target for log(gamma)');
    plot(lg, L);
    xlabel('Log (\gamma)');
    ylabel('L');
end%if DEBUG
end%mne_Jeffreys_prior_estimate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = mne_estimate_scale_noise(Data, Noise, ...
                                                     LeadFields, priorFn)
%MNE_DOUBLE_SEARCH_ESIMATE regularized source estimates with noise scaling
%
% Estimates source covariance for MNE, using a white covariance matrix,
% with a single scaling parameter. 
% Allows an additional parameter to control scaling between noise and data
% as rho*B = LS + e
global DEBUG

sourceCovFn     = @(logGamma) exp(real(logGamma)) .* ...
                              speye(LeadFields.nDims * LeadFields.nSources);
                                    
% put noise and sources are on same scale
noiseToLFScale = median(log(diag(Noise.cov))) - ...
                 median(log(sum(LeadFields.lf.^2)));

% put data and noise on scale set by smallest eigenvalue
dataEigVals      = eig(Data.cov);
noiseEigVals     = eig(Noise.cov);
NoiseToDataScale = log(noiseEigVals(end)) - log(dataEigVals(end));

% set initial values
logGammaInit = noiseToLFScale; 
logRhoInit   = NoiseToDataScale;
paramsInit   = [logGammaInit; logRhoInit];

% priors forms are the same for rho and gamma
priors = [priorFn; priorFn];

% If there is a scale on the priors, we want to constrain the space to
% moderately sensible values
% If we reckon our initial guesses are any good, let's constrain the
% variance parameters to be within a factor of 50=exp(4).
scales = exp(paramsInit + 4);

optimise_target = @(params) optimise_target_double_prior(Data, Noise,        ...
                                                         LeadFields, params, ...
                                                         priors, scales,     ...
                                                         sourceCovFn);

% optimise using a multivariate nonlinear Nelder-Mead minimization
params    = fminsearch(optimise_target, paramsInit);
logGamma  = params(1);
sourceCov = sourceCovFn(logGamma);
rho       = exp(params(2));

% Extract weights for 3d sources
W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf, rho);

if DEBUG,
    lrBound = real(log10(exp(logRhoInit))) + [-7 7];
    lgBound = real(log10(exp(logGammaInit))) + [-15 15];
    nR      = 40;
    nG      = 60;
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
end%if DEBUG
end%mne_Jeffreys_prior_estimate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = sparse_bayes_estimate(Data, Noise, LeadFields, ...
                                                GammaPrior)
%SPARSE_BAYES_ESIMATE regularized source covariance 
%
% Estimates source covariance using sparse Bayes estimate: allow for
% source variance at each location and in each orientation in space. 
%
% To use this algorithm, -log p(gamma) MUST BE CONCAVE for each gamma.
global DEBUG
nSourceElements = LeadFields.nDims * LeadFields.nSources;
sourceCovFn     = @(gamma) spdiags(real(gamma), ...
                                   0, nSourceElements, nSourceElements);
logGammaInit    = (median(log(diag(Noise.cov)))          ...
                   - median(log(sum(LeadFields.lf.^2)))) ...
                  .* ones(nSourceElements, 1);
              
% set scale for priors on gamma
scale = exp(logGammaInit + 5);
              
% use updates from gamma-MAP in Wipf and Nagarajan (2009)
oldGamma   = exp(logGammaInit);
nIterMax   = 1000;
iIter      = 1;
deltaGamma = NaN(nIterMax);

HALT_CONDITION = 1e-5; %fractional change in norm(gamma) to class as convergence

while iIter <= nIterMax,
    if DEBUG && ~mod(iIter, 10),
        fprintf('%s: SB g-MAP iteration no %d.\n', mfilename, iIter);
    end%if
    
    iIter = iIter + 1;
    
    % update rules
    for i = nSourceElements:-1:1,
        sigmaB   = empirical_bayes_cov(Noise.cov, Data.cov, LeadFields.lf);
        gamma(i) = oldGamma(i) ./ sqrt(Data.nSamples) ...
                   *  1            ...
                   * sqrt(trace(LeadFields.lf(:,i).' * inverse(sigmaB) * LeadFields.lf(:,i) ...
                          - GammaPrior.diff_fn(oldGamma(i), scale) ./ Data.nSamples));
    end%for
    
    % monitor change
    deltaGamma(iIter) = norm(gamma - oldGamma)./norm(oldGamma);
    
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

% Extract weights for 3d sources
W = estimate_weights(Noise.cov, sourceCov, LeadFields.lf);

if DEBUG
    figure('Name', 'Sparse-Bayes convergence', 'Color', 'w');
    plot(deltaGamma, 'r', 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', 14);
    ylabel('Fractional change in norm(\gamma)', 'FontSize', 14);
end%if
end%sparse_bayes_estimate





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = optimise_target_single_prior(Data, Noise, LeadFields, ...
                                       logGamma, prior, scale, sourceCovFn)
%OPTIMISE_TARGET_single_PRIOR
% minimise L to find gamma-MAP solution
global DEBUG
persistent iTARGETCALL

Sigma_EB = empirical_bayes_cov(Noise.cov, sourceCovFn(logGamma), ...
                               LeadFields.lf);
                           
% we want to maximise
% log p(g|B) = -0.5 Tr(BB' Sigma_b^{-1}) - n/2 logdet Sigma_b + log p(g).

% use property sum(eig(B, A)) = trace(inv(A) * B)
% or trace(AB) = sum(sum(A .* B')) (and covariance matrices are symmetric)
L = (trace(Data.cov * inverse(Sigma_EB, 'symmetric'))    ...               % faster than elementwise product or sum(eig()). 
     + ROInets.logdet(Sigma_EB, 'chol')) * Data.nSamples ...
    + 2 * prior(logGamma, scale);                           


if DEBUG,
    if ~exist('iTARGETCALL', 'var') || isempty(iTARGETCALL),
        iTARGETCALL = 1;
    else
        iTARGETCALL = iTARGETCALL + 1;
    end%if
    fprintf(['Call to optim fn %4.0d: L = %0.8G, logGamma = %0.4G, ', ...
             'logdet(Sigma_EB) = %0.6G. \n'],                         ...
            iTARGETCALL, L, logGamma, ROInets.logdet(Sigma_EB, 'chol'));
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
pg       = priors(1);
pr       = priors(2);

Sigma_EB = empirical_bayes_cov(Noise.cov, sourceCovFn(logGamma), ...
                               LeadFields.lf);

% we want to maximise
% log p(r,g|B) = -0.5 r^2 Tr(BB' Sigma_b^{-1}) - n/2 logdet Sigma_b  
%                + n*log(r) + log p(r) + log p(g).
%
% use property sum(eig(B, A)) = trace(inv(A) * B)
% or trace(AB) = sum(sum(A .* B')) (and covariance matrices are symmetric)
L = (exp(2*logRho) * trace(Data.cov * inverse(Sigma_EB, 'symmetric')) ...  % faster than elementwise product or sum(eig()). 
     + ROInets.logdet(Sigma_EB, 'chol') - 2*logRho) .* Data.nSamples  ...
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
            ROInets.logdet(Sigma_EB, 'chol'));
end%if DEBUG
end%optimise_target_double_prior





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lp, diff_lp] = parse_priors(priorChoice)
%PARSE_PRIORS creates function handle for priors on hierarchical variance
%parameters. 
%
% Some priors require scales, so lp takes two arguments, lp(loggamma, scale).
% May also require d(lp)/dg. 

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
        
    case 'log-uniform'
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
        Noise.cov = Noise.lamdba .* diag(diagDataCov);
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
%        .normalise       - sLoreta, [norm], none
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
       isfield(x, 'leadFields') && isscalar(leadFields);
IP.addParamValue('ReduceRank', [], test);
test = @(x) ischar(validatestring(x,                                   ...
                                  {'Wens', 'MNE', 'MNE-INA',           ...
                                   'MNE-scaled-noise', 'sparseBayes'}, ...
                                  mfilename, ...
                                  'Options.sourceModel'));
IP.addParamValue('sourceModel', 'MNE', test);
test = @(x) ischar(validatestring(x, {'sLoreta', 'norm', 'none'}, ...
                                  mfilename, ...
                                  'Options.normalise'));
IP.addParamValue('normalise', 'norm', test);
test = @(x) isa(x, 'function_handle') || ...
            ischar(validatestring(x, {'uniform-sd', 'uniform-variance', ...
                                      'log-uniform', 'cauchy',          ...
                                      'normal', 'log-normal'},          ...
                                  mfilename, 'Options.gammaPrior'));
IP.addParamValue('gammaPrior', 'uniform-sd', test);

% parse inputs
IP.parse(Options);

ParsedOptions = IP.Results;

% set prior functions
[ParsedOptions.gamaPrior.fn, ...
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
% Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis 
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
% [EOF]