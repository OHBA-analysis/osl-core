function mats = run_correlation_analysis(nodeData, envData, Regularize)
%RUN_CORRELATION_ANALYSIS runs various correlations on node data
% MATS = RUN_CORRELATION_ANALYSIS(NODEDATA, ENVDATA)
%   produces matrices in structure MATS which give the correlation of the 
%   raw node timecourses, NODEDATA, and the marginal correlation and 
%   partial correlation of the envelope timecourses ENVDATA. 
%
% MATS = RUN_CORRELATION_ANALYSIS(NODEDATA, ENVDATA, REGULARIZE) applies
%   regularization to the estimation of partial correlation matrices using
%   parameters set in REGULARIZE:
%            .do            : true or false: controls use of
%                             regularization. [false]
%            .method        : 'Bayesian' or 'Friedman': use Wang (2012)'s
%                             Bayesian graphical lasso or Friedman (2007)'s
%                             graphical lasso
%            .path          : path of regularization parameters controlling
%                             the strength of the regularization, the best 
%                             value being found using 10-fold CV. 
%                             [Friedman only]
%            .adaptivePath  : true or false: adapt path if the best
%                             regularization parameter is on the edge of
%                             the path [Friedman only] 
%            .Prior         : structure with fields controlling the shape
%                             of the gamma(x; a, 1/b) hyperprior on the 
%                             regularization parameter. 
%                             If this field is not set, the default Kerman 
%                             (2011) neutral hyperprior is used (a=1/3, b=0)
%                             [Bayesian only]
%                .a - shape
%                .b - 1/scale
%                There is no proper uninformative prior which is flat in
%                log x. Common minimally informative priors are a=epsilon,
%                b=epsilon, which is not flat in log x or proper in the
%                limit epsilon -> 0. It can be very informative in datasets
%                with small variances. 
%                The Kerman prior is more uniform in log x and allows the
%                median of the posterior distribution to be specified by
%                the data alone. 
%
% References on techniques used:
%  Friedman, J. and Hastie, T and Tibshirani, R. "Sparse inverse covariance
%  estimation with the graphical lasso", Biostatistics 9(3), 432-441,
%  (2008). 
%  
%  Wang, H. "Bayesian graphical lasso models and efficient posterior
%  computation", Bayesian Analysis 7(2), 771-790 (2012). 
%
%  Kerman, J. "Neutral noninformative and informative conjugate beta and
%  gamma prior distributions", Electronic Journal of Statistics 5,
%  1450-1470 (2011). 
%
%  Gelman, A. "Prior distributions for variance parameters in hierarchical
%  models", Bayesian Analysis 1(3), 515-533 (2006). 


%	Copyright 2014 OHBA
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


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 368 $
%	$LastChangedDate: 2014-12-13 19:05:24 +0000 (Sat, 13 Dec 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 07-Nov-2013 12:43:45

if nargin < 2 || isempty(Regularize), 
    Regularize.do = false; 
end%if 

if ~isempty(nodeData),
    % raw correlations
    rawCorr  = corr(nodeData');
    clear nodeData
else
    rawCorr = [];
end%if

fprintf(' Running correlation analysis on enveloped data: \n');
nodeCov       = cov(real(envData'));
nodeCorr      = corrcov(nodeCov);
nodePrecision = ROInets.cholinv(nodeCov);
nodePCorr     = ROInets.convert_precision_to_pcorr(nodePrecision);

nSamples = ROInets.cols(envData);

mats.nSamples              = nSamples;
mats.correlation           = rawCorr;
mats.envCovariance         = nodeCov;
mats.envCorrelation        = nodeCorr;
mats.envPrecision          = nodePrecision;
mats.envPartialCorrelation = nodePCorr;

if Regularize.do,
    switch Regularize.method 
        case 'Bayesian'
        % do Bayesian Graphical L1 Lasso using Wang (2012)'s method and code
        fprintf('   Regularizing using Bayesian graphical lasso \n');
        % default - use Kerman neutral hyperprior (2011)
        a_lambda = Regularize.Prior.a; % def: 1/3
        b_lambda = Regularize.Prior.b; % def: 0
        
        % MCMC parameters
        burnIn = 3000; %iterations
        nMC    = 8000; %iterations
        
        [~,             ...
         postPrecision, ...
         postLambda] = ROInets.BayesGLasso_Columnwise(nodeCov .* nSamples, ...
                                                      nSamples,            ...
                                                      nodeCov,             ...
                                                      nodePrecision,       ...
                                                      a_lambda,            ...
                                                      b_lambda,            ...
                                                      burnIn,              ...
                                                      nMC);
        
        regPrec  = mean(postPrecision, 3);
        regPCorr = ROInets.convert_precision_to_pcorr(regPrec);
        meanRho  = mean(postLambda) ./ nSamples; % lambda = N * rho
        fprintf('   Mean regularization rho: %g. \n', meanRho);
        
        mats.envPartialCorrelationRegularized = regPCorr;
        mats.envPrecisionRegularized          = regPrec;
        mats.Regularization.mean              = meanRho;
        mats.Regularization.posteriorRho      = postLambda ./nSamples;
        
        
        case 'Friedman'
        % do Friedman (2008)'s graphical lasso. May be faster. 
        fprintf('   Regularizing using graphical lasso and x-validation \n');
        Kfold                  = 10;
        [regPrecision, rhoOpt] = ROInets.glasso_cv(real(normalise_vectors(...
                                           ROInets.demean(envData,2),2)), ...
                                                   Regularize.path,       ...
                                                   Kfold,                 ...
                                                   [],                    ...
                                                   Regularize.adaptivePath);
        
        regPCorr = ROInets.convert_precision_to_pcorr(regPrecision);
        mats.envPrecisionRegularized          = regPrecision;% this is going to be scaled. Sorry.
        mats.envPartialCorrelationRegularized = regPCorr;
        mats.Regularization.mean              = rhoOpt;
        
        otherwise % not one of these methods
        error([mfilename ':UnrecognisedRegMethod'],         ...
              'Unrecognised regularization method: %s. \n', ...
              Regularize.method);
    end%if
    
else
    mats.Regularization.mean = 0;
end%if

end%run_correlation_analysis
                         
 function VV = normalise_vectors(V, dim)
 %NORMALISE_VECTORS normalises rows or columns of a matrix
 %
 % NORMV = GC_NORMALISE_VECTORS(V, DIM) normalises vectors along dimension
 %   DIM of V. If DIM=1, this treats columns as vectors and if DIM=2, this
 %   treats rows as vectors.
 %
 % NORMV = GC_NORMALISE_VECTORS(V) is the same as GC_NORMALISE_VECTORS(V, 1)
 %
 % Example:
 %   If X = [3 3 3]
 %   Then GC_NORMALISE_VECTORS(X,1) is [1 1 1] and GC_NORMALISE_VECTORS(X,2)
 %   is [0.5774 0.5774 0.5774].


 %	Copyright 2013 Giles Colclough
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
 %	$Revision: $
 %	$LastChangedDate: $
 %	Contact: giles.colclough 'at' eng.ox.ac.uk
 %	Originally written on: GLNXA64 by Giles Colclough, 25-Nov-2013 11:49:41

 if nargin < 2 || ~exist('dim', 'var') || isempty(dim),
 dim = 1;
 end%if

 VV = bsxfun(@rdivide, V, sqrt(sum(V.^2, dim) ./ size(V, dim)));
 end%normalise_vectors
% [EOF]
