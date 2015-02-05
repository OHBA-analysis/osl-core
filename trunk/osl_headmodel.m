function D = osl_headmodel(S)
% NOSL_HEADMODEL runs MEG coregistration and forward model in SPM8 or 
% SPM12. These two tasks are separately performed by nosl_datareg.m and 
% osl_forward_model.m which are wrapped together to ensure functionality 
% of local spheres forward model when montaging has been applied
%
% D = osl_headmodel(S)
%
% REQUIRED INPUTS:
%
% S.D               - SPM MEG object filename
%
% S.mri             - structural MRI nii file name (set S.mri=[] or '' to 
%                     use template structural)
%
% S.useheadshape    - set to 0 or 1 to indicated if the headshape points 
%                     should be used in the registration
%
%
% OPTIONAL INPUTS:
%
% S.use_rhino       - use RHINO coregistration instead of SPM
%
% S.forward_meg     - 'Single Shell' or 'MEG Local Spheres' (default)
%
% S.fid_label       - Fiducial labels with fields:
%                       .nasion (Neuromag default 'Nasion')
%                       .lpa    (Neuromag default 'LPA')
%                       .rpa    (Neuromag default 'RPA')
%
% S.fid_mnicoords   - Specify fiducual MNI coordinates with fields:
%                       .nasion - [1 x 3]
%                       .lpa    - [1 x 3]
%                       .rpa    - [1 x 3]
%                       (omit field to use SPM defaults)
% OR:
%
% S.fid_nativecoords - Specify native MNI coordinates with fields:
%                        .nasion - [1 x 3]
%                        .lpa    - [1 x 3]
%                        .rpa    - [1 x 3]
%
% Adam Baker 2014


%%%%%%%%%%%%%%%%%%%%%%%   P A R S E   I N P U T S   %%%%%%%%%%%%%%%%%%%%%%%

% Check SPM File Specification:
try
    S.D = char(S.D);
    [pathstr,filestr] = fileparts(S.D);
    S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
    D = spm_eeg_load(S.D);
catch
    error('SPM file specification not recognised or incorrect');
end%try

% Check Headmodel Specification:
try
    S = ft_checkopt(S,'forward_meg','char',{'Single Shell','MEG Local Spheres'});
catch
    warning('Forward model specification not recognised or incorrect, assigning default: "Single Shell"')
    S = ft_setopt(S,'forward_meg','Single Shell');
end%try

% Check RHINO Specification:
try
    S = ft_checkopt(S,'use_rhino','double');
catch
    S = ft_setopt(S,'use_rhino',0);
end

% Check Structural Specification:
try
    S.mri = char(S.mri);
    [pathstr,filestr,ext] = fileparts(S.mri);   
catch
    S.mri = ft_getopt(S,'mri','');
    error('Structural MRI specification not recognised or incorrect');
end
if ~isempty(S.mri)
    if isempty(ext) % force .nii suffix
        ext = '.nii';
    elseif strcmp(ext,'.gz') && S.use_rhino == 0
        error('S.mri must be .nii (not .nii.gz) when using SPM coregistration')
    else
        tempMesh = spm_eeg_inv_mesh;
        S.mri     = tempMesh.sMRI;
    end
    S.mri = fullfile(pathstr,[filestr,ext]);
end

% Check Headshape Specification:
try
    S = ft_checkopt(S,'useheadshape',{'single','double','logical'},{0,1});
catch
    warning('Headshape specification not recognised or incorrect, assigning default: "1"')
    S = ft_setopt(S, 'useheadshape', 1);
end%try

% Check Fiducial Label Specification:
try
    S = ft_checkopt(S,'fid_label','struct');
    assert(isfield(S.fid_label, 'nasion') &&        ...
           isfield(S.fid_label, 'lpa')    &&        ...
           isfield(S.fid_label, 'rpa'),             ...
           [mfilename ':fid_labelIncorrectFields'], ...
           'Incorrect fields in S.fid_label\n');
catch
    warning('Fiducial label specification not recognised or incorrect, assigning Elekta defaults\n')
    % default
    S = ft_setopt(S,'fid_label',struct('nasion','Nasion', ...
                                       'lpa','LPA',       ...
                                       'rpa','RPA'));
end


if S.use_rhino
    
%%%%%   R U N   C O R E G I S T R A T I O N   U S I N G   R H I N O   %%%%%
    
    S_coreg = S;
    if isfield(S_coreg, 'forward_meg'),
        S_coreg = rmfield(S_coreg, 'forward_meg');
    end%if
    S_coreg.do_plots = 1;
    rhino(S_coreg);
    close all
    
    
%%%%%%%%%%%%%%%%%%   R U N   F O R W A R D   M O D E L   %%%%%%%%%%%%%%%%%%
    
    S_forward               = struct();
    S_forward.D             = S.D;
    S_forward.forward_meg   = S.forward_meg;
    osl_forward_model(S_forward);
    
    
else % ~S.use_rhino
    
%%%%%%%   R U N   C O R E G I S T R A T I O N   U S I N G   S P M   %%%%%%%
    
    matlabbatch{1}.spm.meeg.source.headmodel.D       = {S.D};
    matlabbatch{1}.spm.meeg.source.headmodel.val     = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {[S.mri ',1']};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres    = 2;
    
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = S.fid_label.nasion;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = S.fid_label.lpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = S.fid_label.rpa;
    
    if(isfield(S,'fid_mnicoords'))
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = S.fid_mnicoords.nasion;
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = S.fid_mnicoords.lpa;
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = S.fid_mnicoords.rpa;
    else
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
    end
    
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = S.useheadshape;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = S.forward_meg;
    
    spm_jobman('run', matlabbatch);
    
end % if S.use_rhino 

end%osl_headmodel

