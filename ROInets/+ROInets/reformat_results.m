function correlationMats = reformat_results(mats, Settings)
%REFORMAT_RESULTS	move session correlation mats to frequency band mats
%
% FREQ_MATS = REFORMAT_RESULTS(MATS, SETTINGS) moves session-specific
%   correlation matrices in MATS{iSession}{iFrequency}.correlationMatrix
%   into a new format,
%   FREQ_MATS{iFrequency}.correlationMatrix(:,:,iSession). The SETTINGS
%   structure from oil.ROInetworks must be provided. 

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


%	$LastChangedBy: adambaker86@gmail.com $
%	$Revision: 261 $
%	$LastChangedDate: 2014-10-20 20:19:04 +0100 (Mon, 20 Oct 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 11-Apr-2014 12:14:41

for iSession = Settings.nSessions:-1:1,  
    for iFreq = Settings.nFreqBands:-1:1,
        correlationMats{iFreq}.correlation(:,:,iSession)             = mats{iSession}{iFreq}.correlation;
        correlationMats{iFreq}.envCorrelation(:,:,iSession)          = mats{iSession}{iFreq}.envCorrelation;
        correlationMats{iFreq}.envCovariance(:,:,iSession)           = mats{iSession}{iFreq}.envCovariance;
        correlationMats{iFreq}.envPrecision(:,:,iSession)            = mats{iSession}{iFreq}.envPrecision;
        correlationMats{iFreq}.envPartialCorrelation(:,:,iSession)   = mats{iSession}{iFreq}.envPartialCorrelation;
        correlationMats{iFreq}.envCorrelation_z(:,:,iSession)        = mats{iSession}{iFreq}.env_z;
        correlationMats{iFreq}.envPartialCorrelation_z(:,:,iSession) = mats{iSession}{iFreq}.env_z_partial;
        if isfield(mats{iSession}{iFreq}, 'ARmodel'),
            correlationMats{iFreq}.ARmodel(iSession)                 = mats{iSession}{iFreq}.ARmodel;
        end%if
        if isfield(mats{iSession}{iFreq}, 'H0Sigma'),
            correlationMats{iFreq}.H0Sigma(iSession)                 = mats{iSession}{iFreq}.H0Sigma;
        end%if
        correlationMats{iFreq}.nEnvSamples(iSession)                 = mats{iSession}{iFreq}.nSamples;
        if Settings.Regularize.do,
            correlationMats{iFreq}.envPartialCorrelationRegularized(:,:,iSession)   = mats{iSession}{iFreq}.envPartialCorrelationRegularized;
            correlationMats{iFreq}.envPrecisionRegularized(:,:,iSession)            = mats{iSession}{iFreq}.envPrecisionRegularized;
            correlationMats{iFreq}.envPartialCorrelationRegularized_z(:,:,iSession) = mats{iSession}{iFreq}.env_z_partial_reg;
            correlationMats{iFreq}.Regularization(iSession)                         = mats{iSession}{iFreq}.Regularization;
        end%if
        correlationMats{iFreq}.sessionNames{iSession}                = mats{iSession}{iFreq}.sessionName;
    end%for
end%for

end%reformat_results
