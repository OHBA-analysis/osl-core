function [pcorr, invVar, Aguess] = convert_precision_to_pcorr(precision)
%CONVERT_PRECISION_TO_PCORR converts to partial correlation
%
% PCORR = CONVERT_PRECISION_TO_PCORR(PRECISION) converts PRECISION
%   matrix to partial correlation matrix PCORR. 
%
% [PCORR, INVVAR] = CONVERT_PRECISION_TO_PCORR(PRECISION) outputs the
%   inverse of the variance of each variable, INVVAR. 
%
% See also: CORR, CORRCOV. 

% Additional functionality:
% [PCORR, INVVAR, ESTNETMAT] = CONVERT_PRECISION_TO_PCORR(PRECISION)
%   provides an esimated network matrix under the assumptions of a
%   structural equation model (SEM) z=Az + e for variables z, network 
%   matrix A  and Gaussian noise e. 
%
% See as a reference, Mark W. Woolrich & Klaas E. Stephan. (2013)
%   "Biophysical network models and the human connectome," NeuroImage, 
%   volume 80, pages 330-338, ISSN 1053-8119, 
%   http://dx.doi.org/10.1016/j.neuroimage.2013.03.059.
%  (http://www.sciencedirect.com/science/article/pii/S1053811913003091) 
%  The relvant derivation is in Appendix A. 


%	Copyright 2013 OHBA
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
%	$Revision: 213 $
%	$LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 18-Nov-2013 13:42:28

% find the inverse variances
invVar = diag(precision);

% check for blank input
if isempty(precision),
    pcorr  = [];
    Aguess = [];
    return
end%if

% check for zeros
if any(~invVar),
    warning([mfilename ':InfiniteVariance'], ...
          ['At least one zero was found on the diagonal of the precision ', ...
           'matrix. \nThis is unphysical, corresponding to infinite ', ...
           'variance. There is no corresponding partial correlation ', ...
           'matrix. \n']);
       pcorr = NaN(size(precision));
       Aguess = [];
       return
end%if zeros on diagonal

% construct partial correlation matrix
D = diag(sqrt(1.0./invVar));

L = D * precision * D;
pcorr = - triu(L, 1) - triu(L, 1)' + diag(ones(size(precision,1),1));


% construct esimated SEM network matrix if desired
if nargout>2,
    Aguess = estimate_SEM_network_matrix(precision);
end%if
end%convert_precision_to_pcorr

function estNetMat = estimate_SEM_network_matrix(precision)
% consider an SEM z=Az + e
% (I - A) prop sqrtm(inv(covariance))
estNetMat = - sqrtm(precision); 
% ignore diagonal elements - set to zero
estNetMat = triu(estNetMat, 1) + triu(estNetMat, 1)';
end%estimate_network_matrix
% [EOF]