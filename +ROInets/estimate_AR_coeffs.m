function [AR_coeffs, varEst, AR_partial_coeffs] = estimate_AR_coeffs(data, order)
%ESTIMATE_AR_COEFFS fits an AR model to voxel data and estimates coefficients
% [AR_COEFFS, VAR_EST, AR_PARTIAL_COEFFS] = estimate_AR_coeffs(DATA, P) fits
%   a P-order AR model to every row of DATA, then provides the median
%   coefficients from the fits, AR_COEFFS, the median noise variance VAR_EST
%   and the median AR partial coefficients AR_PARTIAL_COEFFS, useful for
%   model testing.
%
%   AR_COEFFS obey the matlab convention: SUM_{k=0}^{p} AR(k) y(n-k) = e(n)
%   where e(n) is a Gaussian noise process.
%
%	See also nets_r2z, ARYULE

%	References: http://www.math.utah.edu/~zhorvath/ar1.pdf

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
%	Originally written on: GLNXA64 by Giles Colclough, 16-Apr-2014 15:25:19

% Use matlab convention on signs.


fprintf(' Estimating AR model parameters. \n');

% parse input order
if 0 == order, 
    AR_coeffs         = 1;
    varEst            = 0;
    AR_partial_coeffs = 1;
    return
    
elseif order < 0, 
    error([mfilename ':InvalidOrderParameter'], ...
          'The model order must be a non-negative integer. \n');
end%if

nNodes = ROInets.rows(data);

% fit AR model with Yule-Walker method at each row
for iNode = nNodes:-1:1,
    [A(:,iNode), E(iNode), K(:,iNode)] = aryule(data(iNode, :), order);
end%for

% take estimate over all voxels
AR_coeffs         =  median(A,2);
varEst            =  median(E);

% partial autocorrelation can help with model selection pacf = -K
AR_partial_coeffs = -median(K,2);

% check stability of model
if order > 0 && ~is_stable(AR_coeffs),
    error([mfilename ':UnstableModel'], ...
          ['The fitted AR model is unstable. ', ...
           'You could try using more coefficients. \n']);
end%if
end%estimate_AR_coeffs

function isStable = is_stable(AR_coeffs)
% BOOL = IS_STABLE(AR_COEFFS) checks the stability of AR model with
%    coeffecients AR_COEFFS and returns BOOL true or false. 
%
%   The stability of an AR model is confirmed if the roots of the
%   characteristic equation lie within the unit circle. 

% Ref: http://davegiles.blogspot.co.uk/2013/06/when-is-autoregressive-model.html

r        = roots(AR_coeffs);
isStable = all(abs(r) < 1.0);
end%is_stable
% [EOF]