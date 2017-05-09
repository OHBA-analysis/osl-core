%% AFRICA EXAMPLE

% 1. Load in some data

% 2. Run everything manually - make sure 3 components are selected
D_manual = osl_africa(D,'used_maxfilter',1,'artefact_channels',{'EOG','ECG'});

% Results stored in D.ica
D_manual.ica
has_montage(D)
D.n_channels


% If you re-run, existing results are used
D_manual = osl_africa(D_manual,'used_maxfilter',1,'artefact_channels',{'EOG','ECG'});

% Do it automatically - automatically chooses 2 components
D.ica;
D = spm_eeg_load(D.fname);
D_automatic = osl_africa(D,'used_maxfilter',1,'artefact_channels',{'EOG','ECG'},'ident_func',@identify_artefactual_components_auto)

% If you did it automatically, you can edit the selection afterwards
D_touchup = osl_africa(D_automatic,'used_maxfilter',1,'artefact_channels',{'EOG','ECG'});

% Also plotting options for automatic
osl_africa(D2,'used_maxfilter',1,'artefact_channels',{'EOG','ECG'},'ident_func',@identify_artefactual_components_auto,'ident_params',struct('do_mains',true,'do_kurt',true,'do_plots',true));



