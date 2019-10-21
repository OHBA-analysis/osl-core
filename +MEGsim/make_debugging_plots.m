function make_debugging_plots(Dsimulated, ...
                              varNiftiFileName, ...
                              dipLocNiftiFileName, ...
                              beamformedOil, ...
                              ReconResults, ...
                              P, ...
                              brainCoords, ...
                              simulatedDipolePositions, ...
                              NAIniftiFileName, ...
                              simDataReconResults, ...
                              iTrial, ...
                              dipoleSignals, ...
                              LCMVdipoleMag, ...
                              dipMeshInd, ...
                              leadFields, ...
                              doBeamformer)
%MAKE_DEBUGGING_PLOTS generate plots for debugging purposes
%
% re-factored code from main function body of GC_simulate_MEG_data.m

%   Copyright 2014 OHBA
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


%   $LastChangedBy: giles.colclough@gmail.com $
%   $Revision: 213 $
%   $LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%   Contact: giles.colclough 'at' magd.ox.ac.uk
%   Originally written on: GLNXA64 by Giles Colclough, 25-Oct-2013 12:43:45


% Make topoplot of variance in first trial over channels
plot_channel_power(Dsimulated);

% view sensor data
oslview(Dsimulated);

if doBeamformer,
    % view variance in source-space
    fslview({varNiftiFileName; dipLocNiftiFileName}, ...
            [2.3 5; 0 1.2], ...
            {'Red-Yellow'; 'Green'});
        
    % view beamformed variance maps after enveloping
    varMapName = fullfile(beamformedOil.source_recon.dirname, ...
                          beamformedOil.enveloping.name, ...
                          'variance_maps', ...
                          'session1_recon_noise_corrected_variance_map');
    fslview({varMapName; dipLocNiftiFileName}, ...
            [2.3 5; 0 1.2], ...
            {'Red-Yellow'; 'Green'});
end%if

% make plot of brain grid, sensor positions and dipole location to check
% co-ordinate systems
brainMask = nii.quickread(ReconResults.mask_fname, P.spatialRes) > 0;
plot_dipole_positions(brainCoords, ...
                      brainMask, ...
                      sensors(Dsimulated, 'MEG'), ...
                      simulatedDipolePositions, ...
                      ReconResults.BF.data.transforms.toMNI);

% view LCMV result
if doBeamformer,
    fslview({NAIniftiFileName; dipLocNiftiFileName}, ...
            [2.3 5; 0 1.2], ...
            {'Red-Yellow'; 'Green'});
end%if

% Make plots of fit between true signal and reconstructed signal at
% each dipole location
if doBeamformer,
    OILdipoleMag = MEGsim.get_oil_dipole_magnitudes(simDataReconResults, ...
                                                    iTrial);
    % just use first trial
    plot_reconstruction_fits(dipoleSignals{iTrial}.^2, ...
                             LCMVdipoleMag(dipMeshInd,:), ...
                             OILdipoleMag(dipMeshInd,:), ...
                             Dsimulated.time);
end%if

% Make plots of lead fields
% LCMV and simulation use same lead fields.
% oil reconstruction generates a new set, after new source recon stage
if doBeamformer,
    leadFieldPlotLocs = [1 dipMeshInd(1)];
    plot_lead_fields(leadFieldPlotLocs, ...
                     leadFields, ...
                     simDataReconResults.BF.sources.L.MEG);
end%if
end%make_debugging_plots
%%




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_channel_power(Dsimulated)
%PLOT_CHANNEL_POWER plots the power in each channel of D as a topoplot

% calculate power in each channel
iTrial = 1;
channelPower = mean(Dsimulated(:,:,iTrial).^2, 2);

% find the names of all MEG channel types
modalities = MEGsim.find_MEG_chantypes(Dsimulated);

if isempty(modalities),
    error([mfilename ':UnrecognisedMEGchannelTypes'], ...
          ['Neither Neuromag nor CTF channel types appear ', ...
           'to have been simulated. \n']);
end%if

% separate plot for each modality
plotChannels = MEGsim.megchannels(Dsimulated);
for m = 1:length(modalities),
    figure('Color', 'w', 'Name', ...
           sprintf('Power in first trial over sensors for modality %s', ...
                   modalities{m}));
    component_topoplot(Dsimulated, ...
                                 channelPower(plotChannels), ...
                                 modalities{m});
end%loop over modalities
end%plot_channel_power
%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function plot_dipole_positions(MNIcoords, brainMask, grad, dipolePositions, M)
%PLOT_DIPOLE_POSITIONS 
%  plots the masked computation grid, the relative positions of the 
%  MEG sensors, and the location of simulated dipoles. 

internalGrid = MNIcoords(brainMask, :);
sensorGrid   = grad.chanpos;

% we can apply a transformation to MNI
if nargin>4,
internalGrid    = spm_eeg_inv_transform_points(M, internalGrid);
sensorGrid      = spm_eeg_inv_transform_points(M, sensorGrid);
dipolePositions = spm_eeg_inv_transform_points(M, dipolePositions);
end%if

figure('Name', 'dipole grid', 'Color', 'w');
hold on;
pointSize = 4;
% calculation grid
scatter3(internalGrid(:,1), ...
         internalGrid(:,2), ...
         internalGrid(:,3), ...
         pointSize, 'k.');
     
% dipole positions
scatter3(dipolePositions(:,1), ...
         dipolePositions(:,2), ...
         dipolePositions(:,3), ...
         25, 'filled', 'ro');
     
% MEG sensors
scatter3(sensorGrid(:,1), sensorGrid(:,2), sensorGrid(:,3), 'go');

hold off;
xlabel('x');
ylabel('y');
rotate3d;
end%plot_dipole_positions
%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [] = plot_reconstruction_fits(dipoleMag, ...
                                       LCMVdipoleMag, ...
                                       OILdipoleMag, ...
                                       time)
% plots a comparison between simulated and reconstructed signals at the
% dipole locations

nDipoles = size(dipoleMag, 1);

figure('Color', 'w', 'Name', 'OIL beamformer signal fit comparison');
for iDipole = 1:nDipoles,
    subplot(1, nDipoles, iDipole);
    scatter(dipoleMag(iDipole,:), OILdipoleMag(iDipole,:), 'k+');
    set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
    xlabel('dipole signal', 'FontSize', 14);
    if iDipole==1,
        ylabel('Reconstructed signal', 'FontSize', 14);
    end%if
    
    R = corrcoef(dipoleMag(iDipole,:), OILdipoleMag(iDipole,:));
    title(sprintf('Correlation %1.2f', R(1,2)), 'FontSize', 14);
end%loop over dipoles

figure('Color', 'w', 'Name', 'LCMV beamformer signal fit comparison');
for iDipole = 1:nDipoles,
    subplot(1, nDipoles, iDipole);
    scatter(dipoleMag(iDipole,:), LCMVdipoleMag(iDipole,:), 'k+');
    set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
    xlabel('dipole signal', 'FontSize', 14);
    if iDipole==1,
        ylabel('Reconstructed signal', 'FontSize', 14);
    end%if
    
    R = corrcoef(dipoleMag(iDipole,:), LCMVdipoleMag(iDipole,:));
    title(sprintf('Correlation %1.2f', R(1,2)), 'FontSize', 14);
end%loop over dipoles

scale_rows = @(A) bsxfun(@rdivide, A, std(A, 0, 2));

dipoleMagScaled = scale_rows(dipoleMag);
OILscaled       = scale_rows(OILdipoleMag);
LCMVscaled      = scale_rows(LCMVdipoleMag);

% plot actual dipole timeseries
figure('Color', 'w', 'Name', 'Variance-normalised Reconstructed signals');
for iDipole = nDipoles:-1:1,
    hax1(iDipole)  = subplot(2, nDipoles, iDipole);
    set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
    hold on;
    plot(time, dipoleMagScaled(iDipole,:), ...
        'k', 'LineWidth', 2);
    plot(time, OILscaled(iDipole,:), ...
        'Color', [132 112 255]/255, 'LineWidth', 2);
    plot(time, LCMVscaled(iDipole,:), ...
        'Color', [1 0 1], 'LineWidth', 2);
    lsline;
    hold off;
    hLeg = legend('dipole', 'OIL recon', 'LCMV recon');
    set(hLeg, 'FontSize', 12);
    xlabel('Time (s)', 'FontSize', 14);
    axis([0 2 -Inf Inf]);
    if iDipole==1,
        ylabel('Signals', 'FontSize', 14);
    end%if
    
    hax2(iDipole) = subplot(2, nDipoles, iDipole+nDipoles);
    set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
    hold on;
    plot(time, dipoleMagScaled(iDipole,:) - OILscaled(iDipole,:), ...
        'Color', [132 112 255]/255, 'LineWidth', 2);
    plot(time, dipoleMagScaled(iDipole,:) - LCMVscaled(iDipole,:), ...
        'Color', [1 0 1], 'LineWidth', 2);
    lsline;
    hold off;
    xlabel('Time (s)', 'FontSize', 14);
    axis([0 2 -Inf Inf]);
    if iDipole==1,
        ylabel('Residuals', 'FontSize', 14);
    end%if
    
    linkaxes([hax1(iDipole), hax2(iDipole)], 'xy');
end%loop over dipoles
end%plot_reconstruction_fits
%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [] = plot_lead_fields(leadFieldPlotLocs, ...
                               simLeadFields, ...
                               oilLeadFields)
%PLOT_LEAD_FIELDS compare lead fields used for simulation and source
% reconstruction. 

for iLead = leadFieldPlotLocs,
    figure('Name', sprintf('Lead field comparison at voxel %d', ...
                           iLead), ...
           'Color', 'w');
       
    subplot(2,1,1);
    leadsToPlot = simLeadFields{iLead};
    plot(leadsToPlot, 'LineWidth', 2);
    title('Simulation lead field', 'FontSize', 12);
    set(gca, 'FontSize', 14);
    
    subplot(2,1,2);
    leadsToPlot = oilLeadFields{iLead};
    plot(leadsToPlot, 'LineWidth', 2);
    title('OIL recon lead field', 'FontSize', 12);
    set(gca, 'FontSize', 14);
end%for
end%plot_lead_fields
%%
% [EOF]