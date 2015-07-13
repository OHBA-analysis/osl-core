function [D, montage] = osl_montage(S)

% [D, montage] = osl_montage(S)
%
% Wrapper call for spm_eeg_montage that restores the parts of
% D.inv that should still be OK (i.e. it restores the coreg bits, but 
% not the forward model bits). This avoids the need to rerun osl_headmodel,  
% but still requires osl_forwardmodel to be rerun.
%
% MWW & AB 2015

D = spm_eeg_load(S.D);

if isfield(D,'inv')
    D_inv = D.inv;
else
    D_inv=[];
end;

[D,montage] = spm_eeg_montage(S);

if ~isempty(D_inv)
    D.inv = D_inv;
    for ii=1:length(D.inv)
        if isfield(D.inv,'forward')
            D.inv{ii}=rmfield(D.inv{ii},'forward');
        end
        for m = 1:numel(D.inv{ii}.datareg)
            D.inv{ii}.datareg(m).sensors = D.sensors(D.inv{ii}.datareg(m).modality);
        end
    end
end

D.save

