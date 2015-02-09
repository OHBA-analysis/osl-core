function [D, montage] = osl_montage(S)

% [D, montage] = osl_montage(S)
%
% Wrapper call for spm_eeg_montage that restores the parts of
% D.inv that should still be OK (i.e. the coreg bits, but not the
% forward model bits). This avoids the need to rerun osl_headmodel, but
% requires osl_forwardmodel to be rerun.
%
% MWW & AB 2015

D_inv=S.D.inv;

[D, montage] = spm_eeg_montage(S);

D.inv=D_inv;

for ii=1:length(D.inv)
    D.inv{ii}=rmfield(D.inv{ii},'forward');
    D.inv{ii}.datareg.sensors=D.sensors(D.inv{ii}.datareg.modality);
end;

D.save;

