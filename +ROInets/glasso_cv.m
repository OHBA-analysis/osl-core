function [P, rhoBest, W] = glasso_cv(Y, rhoPath, K, criterion, adaptiveGrid)
%GLASSO_CV K-fold cross-validation for shrinkage parameter in glasso
%
% [P, RHOBEST] = glasso_cv(Y, rhoPath, K) performs K-fold
%    cross-validation for the shrinkage parameter rho in a frequentist
%    graphical lasso. RHOBEST is the optimal shrinkage parameter and P the
%    regualised precision matrix using this parameter. Y are nNodes x nSamples
%    data and rhoPath the set of shrinkage parameters to test. 
%	
% [P, RHOBEST] = glasso_cv(Y, rhoPath, K, CRITERION) uses the choice of
%    criterion {'AIC', 'BIC', 'LogLikelihood'} to find the best parameter
%    within the rhoPath. The Bayesian information criterion may tend to
%    select higher regularisation, but there is some work suggesting the
%    corrected Akaike information criterion is the best to use. See
%    references for more. 
%
% [P, RHOBEST] = glasso_cv(Y, rhoPath, K, CRITERION, ADAPTIVEGRID) adapts the
%    rhoPath if RHOBEST is found at the boundary of the previous path when 
%    ADAPTIVEGRID is TRUE. Then the path is finessed three times. 
%
% [P, RHOBEST, W] = glasso_cv(...) also outputs the regularised covariance
%    matrix, W. 
%
%   See also ROInets.glasso_frequentist, ROInets.dp_glasso. 


%   References:
% Foygel "Extended Bayesian Information Criteria for Gaussian
%   Graphical Models"
% Schwarz (1978) "Estimating the dimension of a model"
% Burnham, K. P.; Anderson, D. R. (2004), "Multimodel inference: 
%   understanding AIC and BIC in Model Selection", 
%   Sociological Methods and Research 33: 261?304.
% Friedman (2008) "Sparse inverse covariance estimation with the graphical
%   lasso"
% Menendez (2010) "Gene Regulatory Networks from Multifactorial
%    Perturbations Using Graphical Lasso: Application to the DREAM4 Challenge"
% Mazumder and Hastie, The graphical lasso: New insights and
%   alternatives. Electronic Journal of Statistics 2012. 
% http://en.wikipedia.org/wiki/Akaike_information_criterion
% http://en.wikipedia.org/wiki/Bayesian_information_criterion

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
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 10-Mar-2014 12:09:29

% default settings
if nargin < 3 || isempty(K),
    K = 10;
end%if
if nargin < 4 || isempty(criterion), 
    criterion = 'AIC';
end%if
if nargin < 5 || isempty(adaptiveGrid), 
    adaptiveGrid = true;
end%if

scoreFun = choose_score_fun(criterion);
    
% check valid regression parameters
assert(all(rhoPath >= 0),          ...
       [mfilename, ':InvalidRho'], ...
       'Regression parameters must be non-negative. \n');
assert(all(diff(rhoPath) > 0), ...
       [mfilename, ':DecreasingRhoPath'], ...
       'Regression parameters in rhoPath must be increasing in value. \n');

% check for only one rho value
if 1 == length(rhoPath), 
    S = Y * Y' ./ ROInets.cols(Y); % inline cov
%     P = ROInets.glasso_frequentist(S, rhoPath(1));
    [P, W] = ROInets.dp_glasso(S, [], rhoPath);
    rhoBest = rhoPath;
    return
end%if

% do initial fuction run
[P, rhoBest, W] = main_glasso_cv(Y, rhoPath, K, scoreFun);

% check if we've hit up against the bounds and re-run. Then finesse the
% grid a few times
if adaptiveGrid,
    maxAdaptions = 10;
    rhoMin       = 2 * eps;
    rhoMax       = 150;
    nRho         = length(rhoPath);
    maxFinesses  = 3;
    atPathEdge   = false;
    
    % Change grids if rhoBest is on the boundary of the previous "rhoPath".
    if rhoBest == rhoPath(end),
        % we've hit the upper boundary
        [P, W, rhoBest, rhoPath, atPathEdge] = increase_path(Y, rhoPath,     ...
                                                             rhoBest, K,     ...
                                                             scoreFun, nRho, ...
                                                             rhoMax,         ...
                                                             maxAdaptions,   ...
                                                             atPathEdge);
        
    elseif rhoBest == rhoPath(1),
        % we've hit the lower boundary
        [P, W, rhoBest, rhoPath, atPathEdge] = decrease_path(Y, rhoPath,     ...
                                                             rhoBest, K,     ...
                                                             scoreFun, nRho, ...
                                                             rhoMin,         ...
                                                             maxAdaptions,   ...
                                                             atPathEdge);
    end%if
    
    % finesse the grid
    if ~atPathEdge,
        [P, rhoBest, W] = finesse_path(Y, rhoPath, rhoBest, K, scoreFun, ...
                                       nRho, maxFinesses);
    end%if
end%if

fprintf('      glasso_cv: best rho = %g\n', rhoBest);

end%glasso_cv



%%% SUBFUNCTIONS ----------------------------------------------------------

% scoring functions -------------------------------------------------------
function scoreFun = choose_score_fun(criterion)
%CHOOSE_SCORE_FUN choose scoring criterion
switch criterion
    case 'AIC'
        scoreFun = @akaike_information_criterion;
        
    case 'BIC' 
        scoreFun = @bayesian_information_criterion;
        
    case 'LogLikelihood'
        % take negative as we minimise this function. 
        scoreFun = @(P, S, n) -1 * log_likelihood(P, S, n); 
        
    otherwise

        error([mfilename ':UnrecognisedCriterion'], ...
              'Unrecognised scoring criterion\n');
end%switch
end%choose_score_fun
    
function L = log_likelihood(P, S, n)
% log-likelihood for regularized precision P and covariance S and sample
% size n
if ROInets.isposdef(P),
    L = (n / 2) * (ROInets.logdet(P, 'chol') - trace(S * P));
else
    L = (n / 2) * (ROInets.logdet(P)         - trace(S * P));
end%if
end

function BIC = bayesian_information_criterion(P, S, n)
BIC = -2 * log_likelihood(P, S, n) + log(n) * num_edges(P);
end

function AICc = akaike_information_criterion(P, S, n)
% corrected Akaike information criterion
k    = num_edges(P);
AIC  = -2 * log_likelihood(P, S, n) + 2 * k;
AICc = AIC + 2 * k * (k + 1) / (n - k - 1);
end

function nEdges = num_edges(P, edgeThresh)
% number of edges = model dimension

if nargin < 2 || isempty(edgeThresh),
    % threshold defining a zero in precision matrix
    edgeThresh = 1e-3; % if tolerance of GLASSO algorithm is 1e-4, this seems like a sensible threshold. See Fan et al., 2009
end
uniqueInd = triu(true(size(P)), 1);
nEdges    = sum(abs(P(uniqueInd)) > edgeThresh);
end

% Grid shifting functions -------------------------------------------------
function [P, W, rhoBest, rhoPath, atPathEdge] = increase_path(Y, rhoPath,     ...
                                                              rhoBest, K,     ...
                                                              scoreFun, nRho, ...
                                                              rhoMax,         ...
                                                              maxAdaptions,   ...
                                                              atPathEdge)
% best rho is at top end of grid. 
% Keep increasing grid max until satisfied.

dRho      = rhoPath(end) - rhoPath(end-1); % don't change every iteration otherwise will get huge
nAdaption = 0;

while rhoBest == rhoPath(end),
    fprintf(['      glasso_cv: Moving grids: ', ...
             'current best rho is %g. \n'], rhoBest);

    % check bounds
    if rhoBest >= rhoMax
        warning([mfilename ':MaxRhoValue'], ...
                'Regularization hit max allowed: %g\n', rhoMax);
        atPathEdge = true;
        break
    end%if

    % check iterations
    if nAdaption >= maxAdaptions, 
        warning([mfilename ':MaxPathAdaptions'],                       ...
                'Regularization hit end of path after %d adaptions\n', ...
                maxAdaptions);
        atPathEdge = true;
        break
    end%if

    % move path and re-run
    nAdaption    = nAdaption + 1;
    rhoEnd       = max(nRho * dRho + rhoPath(end), 10 * rhoPath(end));
    rhoPath      = logspace(log10(rhoPath(end-1)), log10(rhoEnd), nRho);

    [P, rhoBest, W] = main_glasso_cv(Y, rhoPath, K, scoreFun);
end%while
end%increase_path

function [P, W, rhoBest, rhoPath, atPathEdge] = decrease_path(Y, rhoPath,     ...
                                                              rhoBest, K,     ...
                                                              scoreFun, nRho, ...
                                                              rhoMin,         ...
                                                              maxAdaptions,   ...
                                                              atPathEdge)
% best rho is at bottom end of grid. 
% Keep decreasing grid max until satisfied.

nAdaption = 0;
while rhoBest == rhoPath(1),
    fprintf('      glasso_cv: Moving grids: current best rho is %g. \n', ...
            rhoBest);
        
    dRho = rhoPath(2) - rhoPath(1); % keep distance in proportion to decreasing rhos
    
    % check bounds
    if rhoBest <= rhoMin, % rhoBest == 0
        rhoBest = max(0, rhoBest);
        atPathEdge = true;
        break
    end%if

    % check iterations
    if nAdaption >= maxAdaptions, 
        warning([mfilename ':MaxPathAdaptions'],                       ...
                'Regularization hit end of path after %d adaptions\n', ...
                maxAdaptions);
        atPathEdge = true;
        break
    end%if

    % move path and re-run
    nAdaption    = nAdaption + 1;            
    rhoStart     = min(rhoPath(1) - dRho, rhoPath(1) / 10);
    if rhoStart < 0, 
        rhoStart = 0; 
    end%if
    rhoPath      = logspace(log10(rhoStart   + eps), ...
                            log10(rhoPath(2)), nRho);

    [P, rhoBest, W] = main_glasso_cv(Y, rhoPath, K, scoreFun);
end%while
end%decrease_path

function [P, rhoBest, W] = finesse_path(Y, rhoPath, rhoBest, K, scoreFun, ...
                                        nRho, maxFinesses)
% test rho on a finer grid between the values on the path surrounding the
% one chosen

nFinesse = 0;
while nFinesse < maxFinesses,
    nFinesse = nFinesse + 1;
    
    fprintf(['      glasso_cv: Finessing grids: level %d \n', ...
             '                  current best rho is %g. \n'], ...
            nFinesse, rhoBest);

    rhoBestInd = find(rhoPath == rhoBest, 1);
    if ~rhoBestInd,
        % we have not found rhoBest in the path. Weird. 
        % keep to the existing value. 
        warning([mfilename ':LostRegParameter'], ...
                ['Regularization parameter mismatch. \n',    ...
                 'Stopping finessing process. \n', ...
                 'Come check me out? \n']);
             
        if 1 == nFinesse, 
            % no results allocated yet
            [P, rhoBest, W] = main_glasso_cv(Y, rhoPath, K, scoreFun);
        end%if
        break
    end%if
    if (1 == rhoBestInd) || (length(rhoPath) == rhoBestInd),
        % we should not be here either, as we had not hit the
        % bounds above. Keep to existing value, again. 
        % unless we've changed which boundary we've hit.
        warning([mfilename ':UnexpectedBoundary'], ...
                ['Regularization parameter unexpectedly hit a ',    ...
                 'path boundary. \nStopping finessing process. \n', ...
                 'Come check me out? \n']);
        if 1 == nFinesse,
            % no results allocated yet
            [P, rhoBest, W] = main_glasso_cv(Y, rhoPath, K, scoreFun);
        end%if
        break
    end%if

    newPathBounds = [rhoPath(rhoBestInd-1), rhoPath(rhoBestInd+1)];
    rhoPath       = logspace(log10(newPathBounds(1)), ...
                             log10(newPathBounds(2)), ...
                             fix(nRho / 1.5));

    [P, rhoBest, W] = main_glasso_cv(Y, rhoPath, K, scoreFun);
end%while
end%finesse_path 

% Main worker function ----------------------------------------------------
function [Pbest, rhoBest, Wbest] = main_glasso_cv(Y, rhoPath, K, scoreFun)
[nNodes, nSamples] = size(Y);

assert(nSamples >= nNodes,                 ...
       [mfilename ':PoorlyFormattedData'], ...
       'Use more measurements than number of nodes\n');
   
k      = floor(nSamples/K);
nTrain = nSamples - k ;
% nTest  = k;

% declare memory
score = zeros(K, length(rhoPath));

% verbosity in loop
verbose = 0;

% run glasso on each fold
ft_progress('init', 'text', '');
for iFold = 1:K
    % create fold
    ft_progress(iFold/K, '      glasso_cv: fold %d out of %d', iFold, K);
    Y_test  = Y(:, ((iFold-1) * k + 1) : (iFold * k));
    Y_train = Y;
    Y_train(:, ((iFold-1) * k + 1) : (iFold * k)) = [];
    
    S_train = Y_train * Y_train' / ROInets.cols(Y_train); % inline cov. Normalise by N not N-1.
    S_test  = Y_test  * Y_test'  / ROInets.cols(Y_test);
    
    P = ROInets.dp_glasso(S_train, [], rhoPath, [], [], verbose);
    for jRho = 1:length(rhoPath)
%         P = ROInets.glasso_frequentist(S_train, rhoPath(jRho), verbose); 
        score(iFold, jRho) = scoreFun(P(:,:,jRho), S_test, nTrain);
    end%for
end%for
ft_progress('close');

% find best fold
[~, rhoBestInd]  = min(mean(score));

% compute output based on best rho, using the whitest of the saved P 
% matrices as a warm start
rhoBest        = rhoPath(rhoBestInd);
S              = Y * Y' / nSamples; % inline cov
[Pbest, Wbest] = ROInets.dp_glasso(S, P(:,:,end), rhoBest);
% Pbest   = ROInets.glasso_frequentist(S, rhoBest);
end%main_glasso_cv
% [EOF]
