function U = cholinv(M, isSparse)
%CHOLINV matrix inverse via cholesky decomposition
%
% U = CHOLINV(M) inverts M using the cholesky decomposition. M must be
%   positive definite. 
%
% U = CHOLINV(M, true) assumes M is sparse and uses the suitesparse set of
%   algorithms. 
%
% If you're going to multiply this matrix onto something else, there are
% better algorithms available. 
%
% If you're looking for robust estimation of an inverse covariance, you
% need to check out the graphical lasso, SCAD, adaptive graphical lasso,
% Bayesian lasso or other Bayesian shrinkage estimator. 


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
%	Originally written on: MACI64 by Giles Colclough, 02-Mar-2015 20:27:04

if 1==nargin || ~isSparse,
    [R,p] = chol(M);
    
    if ~p,
        if 3 == exist('solve_triu', 'file'),
            Rinv = solve_triu(R, eye(size(R,1)));
        else
            Rinv = inv(R); % fast for triu. Faster: solve_triu(R, eye(size(R))); from lightspeed toolbox
        end%if
        U    = Rinv * Rinv';
    else
        warning([mfilename ':NotPosDef'], ...
              'Input matrix not positive definite. Using svd method. \n');
        U = pinv(M);
    end%if
else
    % run sparse algorithms
    U = spinv(M);
end%if
end%cholinv