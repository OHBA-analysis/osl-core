function beta = regression(y, X, nocheck)
%REGRESSION  Solves multivariate regression y = X*b + e using fast mex binaries
%
% BETA = regression(Y, X) solves
%    BETA = pinv(X) * Y <==> y = X * BETA + e, minimising lsq error on e,
%    using a fast qpas mex routine. 
%    If X is not full rank, or Y is a matrix rather than a vector, the 
%    function defaults to Matlab's pinv function
%
% BETA = regression(Y, X, 'nocheck') skips the check on the
%    rank of X, which is quicker if X is known to be full rank in advance. 

%	References:
%	  http://www.stat.colostate.edu/~meyer/hingerev2.pdf
%	

%   Never use beta = X\y: it is much slower than pinv route, and does not
%   handle incomplete rank cases in a sensible manner. 

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
%	$Revision: 214 $
%	$LastChangedDate: 2014-07-24 12:40:42 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 14-Apr-2014 12:24:16


% lock function in memory and define variable which will flag only on first
% run
mlock
persistent PRINT_SPEED_WARNING

if ~exist('PRINT_SPEED_WARNING', 'var') || isempty(PRINT_SPEED_WARNING),
    % This is the first time this function is run this session
    PRINT_SPEED_WARNING = true;
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
    warning([mfilename ':GetMexFiles'], ...
            '% will run much faster using the compiled qpas mex files ',   ...
            'by Adrian Wills. \n',                                         ...
            'They are obtainable under an attribution, non-commercial ',   ...
            'license from \n http://sigpromu.org/quadprog/index.html. \n', ...
            mfilename);
    PRINT_SPEED_WARNING = false; % prevent from displaying on repeated runs
end%if

% check for rank deficiency
if nargin >= 3 && strcmpi(nocheck, 'nocheck'),
    isRankDeficient = false;
else
    isRankDeficient = ROInets.cols(X) > ROInets.rows(X) || rank(X) < ROInets.cols(X);
end%if

% solve problem
if ~useQP || isRankDeficient || ROInets.cols(y) > 1,
    if isRankDeficient,
        fprintf(['%s - X is rank deficient: ',                       ...
                 'using pseudoinverse to calculate regression. \n'], ...
                mfilename);
    end%if
    
    % use pseudoinverse
    beta = pinv(X) * y;
    
else
    % Reformat as quadratic programming problem,
    %    minimise beta' H beta - 2 c' beta
    %    for H = X'X and c = X'y
    
    % note - for y with more than one column, use pinv. For larger
    % problems, that is much faster than looping over columns and using qpas. 
    
    H    =   X' * X;
    f    = - X' * y;
    beta = qpas(H, f);
end%if
end%regression