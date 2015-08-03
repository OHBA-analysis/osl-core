function x = randgamma(shape, scale)
%GC_RANDGAMMA random sample from gamma distribution with shape and scale 
%
% X = GC_RANDGAMMA(A, B) returns a matrix, the same size as A and B, where
%   X(i,j) is a sample from a Gamma(A(i,j), B(i,j)) distribution.
%
%   A and B must be the same size. 
%
%   Gamma(a,b) has density function p(x) = x^(a-1) * exp(-x/b) ...
%                                          / (b^(a) * gamma(a)).
%
%   Mean:     a*b
%   Variance: a*b^2
%   Skewness: 2/sqrt(a)
%   Kurtosis: 6/a
%   Mode:     b*(a-1)
%
%   Same pdf as Matlab's gamrnd
%
%   Tom Minka's lightspeed toolbox must be installed, compiled, and on the
%   Matlab path. 
%
%   See also gamrnd, randg, randgamma. 

%	References:
%	http://en.wikipedia.org/wiki/Gamma_distribution#Scaling


%	Copyright 2014 Giles Colclough
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
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 19-Feb-2014 16:42:31

assert(isequal(size(shape), size(scale)), ...
       [mfilename ':InconsistentInputDimensions'], ...
       'Input dimensions must be the same size. \n');
   
% Check for lightspeed mex file
if 3 == exist('randgamma', 'file'),
    % This code replaces gamrnd(shape, scale) or scale.*randg(shape)
    % using Lightspeed Mex file from a compiled lightspeed toolbox (See Tom
    % Minka's website).
    x = scale .* randgamma(shape);

else
    x = sale .* randg(shape);
end%if

end%GC_randgamma