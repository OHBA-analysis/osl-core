function [Z, z, p] = Fisher_r_to_z(r, N, rho0, sigma)
%FISHER_R_TO_Z	Converts correlations to z-scores
%	
% Z = Fisher_r_to_z(r) uses Fisher's Z-transform to convert
%	correlations r to a Z value. 
%
%   Z-stats will be clipped to a very large value, rather than returning
%   infinity from r==1.
%
% [Z, z, p] = Fisher_r_to_z(r, N) uses Fisher's Z-transform to convert
%	correlations to a Z value, and a z-test statistic for the correlation
%	being different from zero, given a sample size N (Two-tailed, only),
%	assuming a normal null distribution of correlations
%
% [Z, z, p] = Fisher_r_to_z(r, N, rhoTest) compares to an H0:rho =
%   rhoTest instead of H0:rho=0.
%
% [Z, z, p] = Fisher_r_to_z(r, N, rhoTest, SIGMA) uses SIGMA as the
%    standard deviation of the null Z-stat distribution, for example, as
%    found from an empirical (normal) null.

%	References:
%	http://courses.education.illinois.edu/EdPsy580/lectures/correlation-ha.pdf
%	http://support.sas.com/documentation/cdl/en/procstat/63104/HTML/default/viewer.htm#procstat_corr_sect018.htm

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
%	$Revision: 229 $
%	$LastChangedDate: 2014-08-07 20:51:36 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 14-Feb-2014 12:26:45



%%% Constants

% Bound Z to prevent infinities
LARGE_Z = 35; % 35 is a pretty large z-stat!




%%% Input checking

assert(all(abs(r(:))<=1), ...
       [mfilename ':InvalidCorrelations'], ...
       'Correlations must be between -1 and 1. \n');
   
if nargout > 1 && nargin == 1,
    error([mfilename ':NotEnoughInputs'], ...
          'Not enough input arguments\n');
end%if

if nargin > 1, 
    assert(N > 4, ...
           [mfilename ':TooSmallSample'], ...
           'Please use a sample size greater than 4. \n');
end%if

if nargin < 3 || isempty(rho0), 
    rho0 = 0;
else
    assert(abs(rho0) <= 1, ...
           [mfilename ':InvalidTestCorrelation'], ...
           'rho0 must be between -1 and 1. \n');
    if abs(rho0) == 1,
        rho0 = sign(rho0) * tanh(LARGE_Z); % prevent atanh(rho0) being infinity
    end%if
end%if

if nargout > 1 && (nargin < 4 || ~sigma), 
    % default standard deviation under normal assumptions
    sigma = sqrt(1.0 ./ (N-3));
end%if




%%% MAIN     

% Transform correlations to Z
Z = convert_r_to_Z(r, LARGE_Z);

% Do hypothesis testing
if nargout > 1,    
    z = convert_Z_to_z(Z, N, rho0, sigma, LARGE_Z);
    p = ROInets.z_to_p_two_tailed(z); 
end%if
end%Fisher_r_to_z






%%% Subfunctions

function Z = convert_r_to_Z(r, LARGE_Z)
% Fisher's conversion of correlations to Z-distributed statistic
Z           = atanh(r); % correlations of +- 1 will be infinite
Z(isinf(Z)) = LARGE_Z .* sign(Z(isinf(Z)));
end%convert_r_to_Z

function z = convert_Z_to_z(Z, N, rho0, sigma, LARGE_Z)
% construct a normal random variable with mean zero and variance 1/(N-3).
zNoNorm = Z - (repmat(convert_r_to_Z(rho0, LARGE_Z), size(Z)) ...
               - repmat(rho0, size(Z)) ./ (2 * (N-1)));
           
if 0 == sigma, 
    % this can happen if there is very low variance in the orthogonalised
    % data. It leads to the regularisation compressing the empirical
    % distribution to zero. I don't have a good fix for this at the moment.
    sigma = 1.0 ./ sqrt(N-3);
    warning([mfilename ':PoorlyConditionedEmpiricalDistributionWidth'], ...
            ['Width of empirical distribution incorrectly esimated. ',  ...
             'Using normal assumptions. \n']);
end%if

% standard normally distributed variable
z = zNoNorm ./ sigma;
end%convert_Z_to_z
% [EOF]