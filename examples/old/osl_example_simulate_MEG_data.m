%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

global OSLDIR;
    
tilde='/home/mwoolrich/';
tilde='/Users/woolrich/';
osldir=[tilde 'Desktop/osl/osl-core'];    

cd(osldir);
osl_startup()


%% call OSL_EXAMPLE_SIMULATE_MEG_DATA
% example use of MEG data simulation function 

% Giles Colclough 2013

% set dipole positions in MNI coordinates
pos = [- 6,    -30,    60  ; ... 
        18,   -90,   0 ]; 
ori = [];

% generate oscillatory signals for signal definition
Fs       = 100; %Hz
nDipoles = size(pos, 1);
duration = 10; %s
time     = 0:1.0/Fs:duration; % s
 
spatialRes=12; %mm

signalDef = gen_uncorrelated_signals(time, Fs, nDipoles, 'oscillating');

% set template and mri files
templateDir  = '/Users/woolrich/tmpsim/';
templateFile = fullfile(templateDir, 'new250_eo_session1');
%sMRI         = fullfile(templateDir, 'structurals/struct1.nii');

% save simulated object
saveFile = '/Users/woolrich/tmpsim/simulatedMEGdata';

% allow for loading back in the source recon output, for multiple re-runs
% BFresults = load(saveReconFile, 'ReconResultsOut');
% BFresults = BFresults.ReconResultsOut;

% load template
Dtemplate = spm_eeg_load(templateFile);

% load structured noise
% noiseMat = load('/Users/woolrich/tmpsim/neuromag_reduced275_empty_room_noise_covariance');

% simulate data
[Dsimulated, true_dipoleSignals] = osl_simulate_MEG_data(...
                      pos, ...
                      signalDef, ...
                      templateFile, ...
                      'fSample',                      Fs, ...
                      'spatialResolution',            spatialRes, ...  % mm
                      'whiteSignalToNoiseRatio',      0.25, ...
                      'structuredSignalToNoiseRatio', 0, ...  %'emptyRoomNoiseCovariance',     noiseMat.emptyRoomNoiseCovariance, ...
                      'runBeamformer',                true, ...
                      'fileName',                     saveFile, ...
                      'dipoleOrientations',           ori, ...
                      'sourceReconSaveFile',          saveReconFile, ...
                      'modalities',                   'MEGGRAD');


%% Run beamformer on simulated data
  
% Run the source reconstruction (yes, again!), and envelope

% Beamform:    
S                   = [];
S.modalities        = 'MEGGRAD';
S.timespan          = [0 Inf];
%S.pca_order         = 250;
S.type              = 'Scalar';
S.inverse_method    = 'beamform';
S.prefix            = '';
S.dirname           = templateDir;
mni_coords          = osl_mnimask2mnicoords(fullfile(osldir,['../std_masks/MNI152_T1_' num2str(spatialRes) 'mm_brain.nii.gz']));
osl_inverse_model(Dsimulated.fullfile,mni_coords,S);

Dsim_recon=spm_eeg_load(Dsimulated.fullfile);

% save recon to nii file
nii.quicksave(Dsim_recon(:,1:1000,1), fullfile(templateDir, ['recon']));
    
%% correlate recon with dipoleSignals and ouput nii files with those correlations in

clear ccs;
for dd=1:size(dipoleSignals{1},1),
    for ii=1:size(Dsim_recon,1),
        cc=corrcoef(squeeze(Dsim_recon(ii,:,1)),true_dipoleSignals{1}(dd,:));
        ccs(dd,ii)=abs(cc(1,2));
    end
    nii.quicksave(squeeze(ccs(dd,:))', fullfile(templateDir, ['corr' num2str(dd), '_recon']));

end

%% view results

cd(templateDir);
fslview({'recon.nii.gz'; 'corr1_recon.nii.gz';'corr2_recon.nii.gz'});
