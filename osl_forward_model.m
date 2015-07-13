function D = osl_forward_model(S)
% OSL_FORWARD_MODEL runs MEG forward model in SPM8 or SPM12
%
% D = nosl_forward_model(S)
% 
% REQUIRED INPUTS:
%
% S.D             - SPM MEG object filename
%
%
% OPTIONAL INPUTS:
%
% S.forward_meg   - Specify forward model for MEG data 
%                   {'Single Shell' or 'MEG Local Spheres'}
%                     (default 'Single Shell')
% S.forward_eeg   - Specify forward model for EEG data
%
% Adam Baker 2014

% Todo: 
% - Clean up D,Dnew,D2 to make it more obvious what these variables are 





%%%%%%%%%%%%%%%%%%%%%%%   P A R S E   I N P U T S   %%%%%%%%%%%%%%%%%%%%%%%

% Check SPM File Specification:
try
    S.D = char(S.D);
    [pathstr,filestr] = fileparts(S.D);
    S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
    D = spm_eeg_load(S.D);
catch
    error('SPM file specification not recognised or incorrect');
end

% Check Headmodel Specification:
try
    S = ft_checkopt(S,'forward_meg','char',{'Single Shell','MEG Local Spheres'});
catch 
    warning('Forward model specification not recognised or incorrect, assigning default: "Single Shell"')
    S = ft_setopt(S,'forward_meg','Single Shell');
end


% Check coregistration has been run:
if isfield(D,'inv')
    for i = 1:numel(D.inv{1}.datareg)
        if strcmp(D.inv{1}.datareg(i).modality,'MEG')
            D.inv{1}.forward(i).voltype = S.forward_meg;
        elseif strcmp(D.inv{1}.datareg(i).modality,'EEG')
            D.inv{1}.forward(i).voltype = S.forward_eeg;
        end
    end
    D.save;
else 
    error('Coregistration should first be run')
end



%%%%%%%%%%%%%%%%%%   R U N   F O R W A R D   M O D E L   %%%%%%%%%%%%%%%%%%

D = spm_eeg_inv_forward(S.D);
D.save;

end

