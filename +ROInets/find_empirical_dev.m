function [sigma_z,         ...
          sigma_z_partial, ...
          sigma_z_partial_reg] = find_empirical_dev(nIter,           ...
                                                    nNodes,          ...
                                                    nSamples,        ...
                                                    ARmodel,         ...
                                                    doRegularize,    ...
                                                    rho,          ...
                                                    Fs,              ...
                                                    FilterSettings,  ...
                                                    EnvelopeParams)
%FIND_EMPIRICAL_DEV  Build up correlations from empirical null distribtion


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
%	$Revision: 366 $
%	$LastChangedDate: 2014-12-13 19:03:44 +0000 (Sat, 13 Dec 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 16-Apr-2014 17:09:11

if nIter < 1, 
    error([mfilename ':BadnH0Iter'], ...
          'Specify how many simulated data sets to create in nH0Iter');
end%if


if doRegularize,
    assert(all(rho(:) >= 0),                       ...
           [mfilename ':BadRegularizationParameter'], ...
           'Regularization parameters must be greater than zero. \n');
    
    % check conditioning of rho
    if isscalar(rho), 
        rho = repmat(rho, 1, nIter);
    else
        assert(isequal(length(rho), nIter),                      ...
               [mfilename ':NumberRegularizationParameters'],       ...
               ['Please use as many regularization parameters as ', ...
                'iterations, or just set one. \n']);
    end%if
    
else
    rho = 0;
end%if


% produce AR random data
randData = filter(1,              ...
                  ARmodel.coeffs, ...
                  sqrt(ARmodel.varianceEstimate) .* randn(nSamples, nNodes, nIter));
t        = (0:1:nSamples-1) ./ Fs;

for iIter = nIter:-1:1, 
    % bandpass filter the data
    if isempty(FilterSettings.band),
        filteredData = transpose(randData(:,:,iIter));
    else
        filteredData = ft_preproc_bandpassfilter(transpose(randData(:,:,iIter)), ...
                                                 Fs,                             ...
                                                 FilterSettings.band,            ...
                                                 FilterSettings.order,           ...
                                                 FilterSettings.type,            ...
                                                 FilterSettings.direction);
    end%if
    
    % envelope the data
    testData = ROInets.envelope_data(filteredData, t, EnvelopeParams)';
    clear filteredData
    
    % take covariance
    rCov  = real(cov(testData));
    rCorr = corrcov(rCov,1);
    
    % we've had a weird error with infinities. catch that
    isOutOfBounds = any(isinf(rCov(:))     | isnan(rCov(:)))    || ...
                    any(isinf(testData(:)) | isnan(testData(:)));
    if isOutOfBounds, 
        error([mfilename 'CorrelationOutOfBounds'],          ...
              ['Something''s wrong here. check me out. \n ', ...
               'Possibly application of AR model to random data is failing. \n']); 
    elseif ~ROInets.isposdef(rCov),
        c = rCov - rCov';
        if max(abs(c(:))) > eps,
            error([mfilename ':AsymmetricCorrelation'], ...
                  'We have generated an asymmetric covariance matrix. \n');
        else
            warning([mfilename ':NonPosDefCorrelation'], ...
                    'We have generated a non-pos-def covariance matrix. \n');
        end%if
    end%if
    
    clear testData
    
    % extract only unique values
    uniqueInd           = triu(true(size(rCorr)), 1);
    rEmpirical(:,iIter) = rCorr(uniqueInd);
    
    % repeat for partial correlations
    r_partial                   = ROInets.convert_precision_to_pcorr(pinv(rCov));
    rEmpirical_partial(:,iIter) = r_partial(uniqueInd);
    
    % repeat for regularized partial correlations
    if doRegularize
        rPrecision = ROInets.dp_glasso(rCov, [], rho(iIter)); 
        r_reg      = ROInets.convert_precision_to_pcorr(rPrecision);
        
        rEmpirical_partial_reg(:,iIter) = r_reg(uniqueInd);
    else
        rEmpirical_partial_reg          = [];
    end
end%for

sigma_z             = std(ROInets.Fisher_r_to_z(rEmpirical(:)));
sigma_z_partial     = std(ROInets.Fisher_r_to_z(rEmpirical_partial(:)));
sigma_z_partial_reg = std(ROInets.Fisher_r_to_z(rEmpirical_partial_reg(:)));

end%find_empirical_dev
% [EOF]