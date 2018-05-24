function signals = gen_uncorrelated_signals(time, Fs, nDipoles, method)
% wraps up several different methods for generating uncorrelated signals. 

% generate dipole signals
switch lower(method)
    case 'white'
        % simulate some white noise with large power at signal voxels.
        variance = 1;
        nSamples = numel(time);
        for iDipole = nDipoles:-1:1,
            rng('shuffle');
            signals(iDipole, :) = sqrt(variance) .* randn(1,nSamples);
        end% loop over dipoles
        
        % bandfilter
        % zero phase iir filter to extract band of interest
        fBand = [1 60]; % Hz
        filterOrder = 5;
        filterType = 'but'; % butterworth
        signals = ft_preproc_bandpassfilter(signals, ...
                                            Fs, ...
                                            fBand, ...
                                            filterOrder, ...
                                            filterType);        
    case 'oscillating'
        assert(nDipoles<=2, ...
               'Not done for more that 2 dips yet!\n');
        freqs = [0.2; 1];
        phase = [0; pi/4];
        for iDip = nDipoles:-1:1,
            signals(iDip, :) = cos(2*pi*freqs(iDip)*time + phase(iDip));
        end
    otherwise
        error('Method not recognised! \n');
end%switch
% demean
signals = bsxfun(@minus, signals, mean(signals, 2));

% orthogonalise
signals = GC_symmetric_orthogonalise(signals')';

signals = bsxfun(@minus, signals, mean(signals, 2)); % again!
end

function L = GC_symmetric_orthogonalise(A)
%GC_SYMMETRIC_ORTHOGONALISE closest orthogonal matrix
% 
% L = GC_SYMMETRIC_ORTHOGONALISE(A) returns orthogonal matrix L which
%   is closest to A, as measured by the Frobenius norm of (L-A). 
%
%   The orthogonal matrix is constructed from a singular value decomposition
%   of A. 
%
%   See also: GC_MGS_ORTHOGONALISE, SVD. 


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


%	$LastChangedBy: GilesColclough $
%	$Revision: 383 $
%	$LastChangedDate: 2014-03-17 15:26:51 +0000 (Mon, 17 Mar 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45


[U, ~, V] = svd(A, 'econ');
L         = U * conj(transpose(V));
end%GC_symmetric_orthogonalise