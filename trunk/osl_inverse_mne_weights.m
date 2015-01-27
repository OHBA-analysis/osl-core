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
%        .sourceModel    - Wens, MNE, sparseBayes
%        .normalise - sLoreta, norm, none
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

%%%%% Questions for mark:
%%%%%                   why demean the covariance matrices?
%%%%%                   why set the default value of lambda to 1.0./mean(diag(SensorData.cov))?
%%%%%                   why estimate source directions from the inverse
%%%%%                   covariance?
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


%% Parse noise matrix
% noise covariance can be specified as the identity, the diagonal of the 
% data covariance, or passed in as a measured or estimated quantity. 
% For the first two cases, a global scaling is possible with parameter
% lambda. 
Noise = parse_noise(Noise, diag(SensorData.cov));
    
%% Source model
switch lower(Options.sourceModel)
    case 'wens'
        sourceCov = wens_estimate(SensorData, Noise, LeadFields);
    case 'mne'
        sourceCov = mne_estimate(SensorData, Noise, LeadFields);
    case 'sparsebayes'
        sourceCov = sparse_bayes_estimate(SensorData, Noise, LeadFields);
    otherwise
        error([mfilename ':InvalidSourceMethod'], ...
              'The chosen source method %s is not recognised. \n', ...
              Options.method);
end%switch

%% Extract weights using Tikhonov regularised form
% equation 13 from Wipf and Nagarajan (2009).
W_3d = sourceCov * LeadFields.lf.' / ...               % A / B = A * inv(B)
       empirical_bayes_cov(Noise.cov, sourceCov, LeadFields.lf);
   
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

end%empirical_bayes_cov





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sourceCov = wens_estimate(SensorData, Noise, LeadFields)
%WENS_ESTIMATE regularized source covariance by Vincent Wens' method
% regularization parameter from Wens et al. Sec 2.4

% use property sum(eig(B, A)) = trace(inv(A) * B)
k = sum(eig(LeadFields.lf * LeadFields.lf.', Noise.cov)) ./ ...
     (sum(eig(SensorData.cov, Noise.cov)) - SensorData.nSensors); 

% In our formulation, gamma = 1/k and C = I
sourceCov = (1.0 ./ k) * eye(LeadFields.nDims * LeadFields.nSources);

end%wens_estimate





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sourceCov = mne_estimate(Data, Noise, LeadFields)
%MNE_ESIMATE regularized source covariance 
%
% Estimates source covariance for MNE, using a white covariance matrix,
% with a single scaling parameter. Places a Jeffrey's prior on the scale
% parameter, i.e. p(gamma) ~ 1/gamma, or p(log(gamma)) ~ 1.
% Eq 7 from Wipf and Nagarajan (2009)
% 1/gamma = 0.5 exp(-f(gamma)) => f(gamma) = log(gamma) - log(2)
global DEBUG

sourceCovFn     = @(logGamma) exp(logGamma) .* ...
                              eye(LeadFields.nDims * LeadFields.nSources);
                          
optimise_target = @(logGamma) optimise_target_Jeffreys_prior(Data, Noise, ...
                                        LeadFields, logGamma, sourceCovFn);
                                    
% variance of data can change by factor of e^5 relative to noise
noiseToLFScale = median(log(diag(Noise.cov))) - ...
                 median(0.5.*log(sum(LeadFields.lf.^2)));
logGammaBound  =noiseToLFScale + [-10 20]; 

% optimise using a golden section search and parabolic interpolation
logGamma      = fminbnd(optimise_target, logGammaBound(1), logGammaBound(2));
sourceCov     = sourceCovFn(logGamma);

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
function sourceCov = sparse_bayes_estimate(Data, Noise, LeadFields)
%MNE_JEFFREYS_PRIOR_ESIMATE regularized source covariance 
%
% Estimates source covariance using sparse Bayes estimate: allow for
% source variance at each location and in each orientation in space, but
% place a 1/gamma ard prior on each variance parameter
% i.e. p(gamma) ~ 1/gamma, or p(log(gamma)) ~ 1.
% Eq 7 from Wipf and Nagarajan (2009)
%  f(gamma) = log(gamma) - log(2)

sourceCovFn     = @(logGamma) diag(exp(logGamma));

optimise_target = @(logGamma) optimise_target_Jeffreys_prior(Data, Noise, ...
                                        LeadFields, logGamma, sourceCovFn);
                                    
logGammaInit    = (median(log(diag(Noise.cov)))               ...
                   - median(0.5.*log(sum(LeadFields.lf.^2)))) ...
                  .* ones(LeadFields.nDims * LeadFields.nSources, 1);
              
% optimise using a multivariate nonlinear Nelder-Mead minimization
logGamma      = fminsearch(optimise_target, logGammaInit);
sourceCov     = sourceCovFn(logGamma);
end%mne_Jeffreys_prior_estimate





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = optimise_target_Jeffreys_prior(Data, Noise, LeadFields, ...
                                            logGamma, sourceCovFn)
%OPTIMISE_TARGET_JEFFREYS_PRIOR
% minimise L to find gamma-MAP solution
global DEBUG
persistent iTARGETCALL

Sigma_EB = empirical_bayes_cov(Noise.cov, sourceCovFn(logGamma), ...
                               LeadFields.lf);

% use property sum(eig(B, A)) = trace(inv(A) * B)
% or trace(AB) = sum(sum(A .* B')) (and covariance matrices are symmetric)
L = trace(Data.cov * inverse(Sigma_EB, 'symmetric')) ...
    + ROInets.logdet(Sigma_EB, 'chol')               ...
    + (sum(logGamma) - log(2)) ./ Data.nSamples;                           % faster than elementwise product or sum(eig()). 


if DEBUG,
    if ~exist('iTARGETCALL', 'var') || isempty(iTARGETCALL),
        iTARGETCALL = 1;
    else
        iTARGETCALL = iTARGETCALL + 1;
    end%if
    fprintf('Call to optim fn %4.0d: L = %0.8G, logGamma = %0.4G, logdet(Sigma_EB) = %0.6G. \n', ...
            iTARGETCALL, L, logGamma, ROInets.logdet(Sigma_EB, 'chol'));
end%if DEBUG
end%optimise_target_Jeffreys_prior





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

nSensors = length(diagDataCov);

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
% [EOF]