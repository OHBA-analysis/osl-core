function centredX = demean(X, dim)
%DEMEAN	Remove mean value
% CENTRED_X = DEMEAN(X) removes the column means from X or works along
%   first non-zero dimension
%
% CENTRED_X = DEMEAN(X, DIM) removes means from X along dimension DIM.
%	
%	See also MEAN, BSXFUN. 


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
%	Originally written on: GLNXA64 by Giles Colclough, 11-Dec-2013 15:31:14

if nargin == 1,
    if(ROInets.rows(X) > 1)
      dim = 1;
    elseif(ROInets.cols(X) > 1)
      dim = 2;
    elseif isscalar(X)
        centredX = 0; 
        return
    else
        dim = find(size(X), 1, 'first');
    end%if
end%if dim not provided

centredX = bsxfun(@minus, X, sum(X,dim) ./ size(X,dim)); % inline mean(X, dim);

% To test speeds:
%bsxfunElapsed time is 16.242885 seconds.
% repmatElapsed time is 22.516901 seconds.
% y = magic(10000);
% dim = 2;
% 
% fprintf('bsxfun')
% tic;
% for i = 1:100,
%     z = bsxfun(@minus, y, mean(y,dim));
% end
% toc
% 
% fprintf('repmat')
% tic;
% for i = 1:100,
%     dims = size(y);
%     dimsize = size(y,dim);
%     dimrep = ones(1,length(dims));
%     dimrep(dim) = dimsize;
%     z = y - repmat(mean(y,dim), dimrep);
% end
% toc
% 
% 
% end