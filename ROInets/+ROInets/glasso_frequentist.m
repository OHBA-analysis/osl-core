function [P, W] = glasso_frequentist(S, rho, verbose)
%GLASSO_FREQUENTIST	graphical lasso for regularized precision matrix estimation
%	
% P = GLASSO_FREQUENTIST(S, RHO) outputs the regularized precision
%   matrix, P, from the sample covariance matrix S, using L1 norm
%   regularization parameter RHO. 
%
% [P, W] = GLASSO_FREQUENTIST(S, RHO) also returns the regularised
%   covariance matrix W. 
%
% P = GLASSO_FREQUENTIST(S, RHO, VERBOSE) switches between three levels
%   of verbosity: choose between 0,1,2. 
%	See also glasso_FTH, L1precisionBCD, glasso. 

%	References:
%	Jerome Friedman, Trevor Hastie and Robert Tibshirani Sparse inverse covariance estimation with the graphical lasso. Biostatistics, December 12, 2007
%
%   http://www.di.ens.fr/~mschmidt/Software/L1precision.html
%   http://statweb.stanford.edu/~tibs/glasso/index.html
%	http://www.stat.sc.edu/~wang345/RESEARCH/Bglasso/bglasso.html

%   This software builds heavily on contributions by Mark Schmidt
%   (L1precisionBCD) and Hao Wang (glasso_FTH). See references above.
%   Neither code package was released with a license, beyond author
%   acknowledgment. 

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


% QPAS licence (called as subfunction):
%     +----------------------------------------------+
%     | Written by Adrian Wills,                     |
%     |            School of Elec. Eng. & Comp. Sci. |
%     |            University of Newcastle,          |
%     |            Callaghan, NSW, 2308, AUSTRALIA   |
%     |                                              |
%     | Last Revised  25 May 2007.                   |
%     |                                              |
%     | Copyright (C) Adrian Wills.                  |
%     +----------------------------------------------+
%    
%   The current version of this software is free of charge and 
%   openly distributed, BUT PLEASE NOTE:
%   
%   This software must be referenced when used in a published work.
%   
%   This software may not be re-distributed as a part of a commercial product. 
%   If you distribute it in a non-commercial products, please contact me first, 
%   to make sure you ship the most recent version.
%   
%   This software is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%   
%   IF IT FAILS TO WORK, IT'S YOUR LOSS AND YOUR PROBLEM.

%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 214 $
%	$LastChangedDate: 2014-07-24 12:40:42 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 19-Feb-2014 16:56:06

% lock function in memory and define variable which will flag only on first
% run
mlock
persistent PRINT_SPEED_WARNING

if ~exist('PRINT_SPEED_WARNING', 'var') || isempty(PRINT_SPEED_WARNING),
    % This is the first time this function is run this session
    PRINT_SPEED_WARNING = true;
end%if


% check if no regularisation
if 0 == rho,
    if nargout > 1,
        W = S;
    end%if
    P = pinv(S);
    return
elseif ~isscalar(rho),
    error([mfilename ':NonScalarRho'],                        ...
          '%s: Expected scalar regularization parameter. \n', ...
          mfilename); 
end%if

% Check for qp mex file
if exist('qpas', 'file') == 3,
    useQP = 1;
    PRINT_SPEED_WARNING = false;
else 
    useQP = 0;
end%if

% Tell the user they should get compiled mex files
if PRINT_SPEED_WARNING,
    warning([mfilename ':GetMexFiles'],                                    ...
            ['%s will run much faster using the compiled qpas mex files ', ...
             'by Adrian Wills. \n',                                        ...
             'They are obtainable under an attribution, non-commercial ',  ...
             'license from \n ',                                           ...
             '<a href="http://sigpromu.org/quadprog/index.html">',         ...
             'http://sigpromu.org/quadprog/index.html</a>. \n'],           ...
            mfilename);
    PRINT_SPEED_WARNING = false; % prevent from displaying on repeated runs
end%if

% switch off display
if nargin < 3, verbose = 1; end

if verbose == 2, 
    fprintf('    Running glasso:\n');
end%if

optTol = 1e-5;

if useQP % use method from L1precisionBCD, a little stripped down
    
    p       = size(S,1);
    maxIter = 10;
    A       = [eye(p-1,p-1); -eye(p-1,p-1)];
    f       = zeros(p-1,1);
    
    % Initial W
    W = S + rho * eye(p,p);
    
    % qp parameters
    qpSolver = @qpas;
    qpArgs   = {[],[],[],[]};
    
    for iter = 1 : (maxIter+1),
        % Check Primal-Dual gap
        P   = ROInets.cholinv(W); %W ^ -1; % W should be PD
        gap = trace(S*P) + rho * sum(sum(abs(P))) - p;
        if verbose == 2,
            fprintf('      Iter = %d, OptCond = %.5f\n', iter, gap);
        end%if
        if gap < optTol
            if verbose,
                fprintf('   glasso: Solution Found\n');
            end%if
            return
        end%if
        
        for i = 1:p
            if verbose == 2, % this has the potential to cut short a run over all columns ... not necessarily a good idea...
                P   = pinv(W); %W^-1; % W should be PD
                gap = trace(S*P) + rho * sum(sum(abs(P))) - p;
                fprintf('        Column = %d, OptCond = %.5f\n', i, gap);
                if gap < optTol
                    fprintf('   glasso: Solution Found\n');
                    break
                end%if
            end%if
            
            % Compute Needed Partitions of W and S
            s_12 = S(ROInets.setdiff_pos_int(1:p, i), i);
            
            H = 2 * pinv(W(ROInets.setdiff_pos_int(1:p,i), ROInets.setdiff_pos_int(1:p,i)));
            b = rho * ones(2 * (p-1), 1) + [s_12; -s_12];
            w = qpSolver((H + H')/2, f, A, b, qpArgs{:});
            
            % Un-Permute
            W(ROInets.setdiff_pos_int(1:p,i), i) = w;
            W(i, ROInets.setdiff_pos_int(1:p,i)) = w';
        end%for
    end%for
    
else % use Wang's matlab method, as it seems to be a little faster
    
    % glasso_FTH(sample covariance matrix, rho: L1 regulariztion parameter)
    % called as: P = glasso_FTH(S, rho);
    
    % N.B.  This algorithm may collapse in certain cases, especially when p large, as no positive
    % definite contraints are imposed in the glasso algorithm.
    
    % Written by Hao Wang @ U of South Carolina
    
    p = size(S,1);
    if isscalar(rho),
        rho = rho*ones(p);
    end%if
    
    avgSoff = mean(S(~eye(p))); % average of off-diagonal elements of empirical covariance matrix S
    
    t = 1e-4; % a fixed threshold ;
    
    % Initial value
    W = S + diag(diag(rho));
    
    % Maximum number of iterations
    Max1 = 30; % across column
    Max2 = 30; % Within column, gradient descend
    
    for iter1 = 1:Max1
        if verbose == 2,
            fprintf('      Iter = %d\n',iter1);
        end%if
        
        W_old = W;
        
        for i = 1:p
            if i == 1
                ind_noi = (2:p)';
            elseif i==p
                ind_noi = (1:p-1)';
            else
                ind_noi = [1:i-1, i+1:p]';
            end%if
            
            V      = W(ind_noi,ind_noi);
            s12    = S(:,i); 
            s12(i) = [];
            w12    = W(ind_noi,i);
            
            beta = V \ w12;
            
            % below Pathwise coordinate descent
            for iter2 = 1:Max2
                if verbose == 2,
                    fprintf('        Column = %d, OptTol = %.5f\n', ...
                            iter2, optTol);
                end%if
                
                beta_old = beta;
                % Speed-up trial which didn't work!
                %        %%% GC replaces code below
                %        j = (1:p-1)';
                %        x = s12(j) - V(j,:)*beta + diag(V(j,j)) .* beta(j);
                %
                %        signx = sign(x);
                %        signx(x==0) = rand(sum(x==0),1); %- this line was being called lots. Led to massive time increase and warnings about poor scaling in L60 beta = V\w12.
                %
                %        beta(j) = max(0, abs(x) - rho(i, ind_noi(j))') .* signx ./ diag(V(j,j));
                %        %---replaces this 
                for j = 1:p-1
                    x = s12(j) - V(j,:) * beta + V(j,j) * beta(j);
                    if 0 == x,
                        signx = rand(1);
                    else
                        signx = sign(x);
                    end%if
                    
                    beta(j) = max(0, abs(x) - rho(i,ind_noi(j))) ...
                             * signx / V(j,j);
                end%for
                % END of GC trial
                
                if max(abs(beta_old - beta)) < optTol,
                    break
                end%if
            end%for
            
            w12          = V * beta;
            W(i,ind_noi) = w12';
            W(ind_noi,i) = w12;
        end%for
        
        chg = abs(W_old - W);
        
        if mean(chg(:)) < (t * avgSoff),
            break
        end%if
    end%for
    
    P = ROInets.cholinv(W);
    
    if verbose,
        fprintf('   glasso: Solution Found\n');
    end%if
end%if
end%glasso_frequentist
% [EOF]