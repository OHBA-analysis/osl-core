function [ReconResults,dipoleOrientations] = compute_recon(D_orig,structural_files)
    % Take in a D object, where coregistration has already been performed
    original_dir = pwd;

    tempfile = tempname('./_temp_beamforming');
    ra.meg.copy(D_orig.fullfile,tempfile)
    D = spm_eeg_load(tempfile)


    % Update the structurals
    f = fields(D.inv{1}.mesh);
    for j = 1:length(f)
        if ischar(D.inv{1}.mesh.(f{j})) && strfind(D.inv{1}.mesh.(f{j}),'/home/disk3/abaker/Data/Schizophrenia/structurals/')
            D.inv{1}.mesh.(f{j}) = strrep(D.inv{1}.mesh.(f{j}),'/home/disk3/abaker/Data/Schizophrenia/structurals/','/Users/romesh/oxford_postdoc/data/adam/Schizophrenia/structurals/');
        end
    end
    % Change the montage
    D = D.montage('switch',0);
    D.save

    if nargin < 2 || isempty(structural_files) 
        structural_files = '/Users/romesh/oxford_postdoc/data/adam/Schizophrenia/structurals/1001_MPRAGE_flat';
    end

    spm_jobman('initcfg');
    startup.osl

    oat=[];
    oat.source_recon.dirname = D.fullfile;
    oat.source_recon.D_continuous={D.fullfile};
    oat.source_recon.conditions={'Undefined'};
    oat.source_recon.freq_range=[4 30]; % frequency range in Hz
    oat.source_recon.time_range=[0 Inf];
    oat.source_recon.modalities = 'MEG';

    % beamformer parameters.
    oat.source_recon.method='beamform';
    oat.source_recon.gridstep=8; % in mm, using a lower resolution here than you would normally, for computational speed
    oat.source_recon.mri=structural_files;

    oat.to_do=[1 0 0 0];
    oat.source_recon.pca_dim=270;
    % oat.source_recon.dirname=['tempdir.oat'];

    try
    oat = osl_run_oat(oat);
    end

    a = load('session1_recon');
    cd(original_dir)
    ReconResults.BF = a.oat_stage_results.BF
    global eta
    dipoleOrientations = eta;

