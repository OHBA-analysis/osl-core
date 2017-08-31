function D = osl_forward_model(D,varargin)
    % OSL_FORWARD_MODEL runs MEG forward model in SPM8 or SPM12
    %
    % D = osl_forward_model(D,S)
    % 
    % REQUIRED INPUTS:
    % - D - MEEG object
    %
    % OPTIONAL INPUTS:
    % - forward_meg - Specify forward model for MEG data {'Single Shell' or
    %   'MEG Local Spheres'} (default 'Single Shell')
    % - forward_eeg - Specify forward model for EEG data
    %
    % Returns an MEEG object with D.inv contents updated. NOT saved to disk
    %
    % Romesh Abeysuriya 2017 
    % Adam Baker 2014

    arg = inputParser;
    arg.KeepUnmatched = true;
    arg.addParameter('forward_meg','Single Shell',@(x) any(strcmp(x,{'Single Shell','MEG Local Spheres'})));
    arg.addParameter('forward_eeg','EEG BEM');
    arg.parse(varargin{:});

    if ~isfield(D,'inv')
        error('Coregistration should first be run')
    end

    for i = 1:numel(D.inv{1}.datareg)
        if strcmp(D.inv{1}.datareg(i).modality,'MEG')
            D.inv{1}.forward(i).voltype = arg.Results.forward_meg;
        elseif strcmp(D.inv{1}.datareg(i).modality,'EEG')
            D.inv{1}.forward(i).voltype = arg.Results.forward_eeg;
        end
    end

    D = spm_eeg_inv_forward(D);



