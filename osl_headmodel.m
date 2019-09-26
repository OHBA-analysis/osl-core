function D = osl_headmodel(S)
% OSL_HEADMODEL runs MEG coregistration and forward model in SPM8 or 
% SPM12. These two tasks are separately performed by nosl_datareg.m and 
% osl_forward_model.m which are wrapped together to ensure functionality 
% of local spheres forward model when montaging has been applied
%
% D = osl_headmodel(S)
%
% REQUIRED INPUTS:
%
% S.D               - SPM MEG object filename or MEEG
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
% S.forward_eeg     - 'EEG BEM' (default)
%
% S.fid              - Fiducial definition: [] for manual placement, or
%                      define coordinates using the following fields:
%
%                      .label    - Fiducial labels with fields:
%                                   .nasion (Neuromag default 'Nasion')
%                                   .lpa    (Neuromag default 'LPA')
%                                   .rpa    (Neuromag default 'RPA')
%
%                      .coords   - Specify fiducual coordinates with fields:
%                                   .nasion - [1 x 3]
%                                   .lpa    - [1 x 3]
%                                   .rpa    - [1 x 3]
%                                   (leave empty to use SPM defaults)
%
%                      .coordsys - Specify fiducial coordinate system as:
%                                  'Native' or 'MNI' (default 'MNI')
%
% Romesh Abeysuriya 2017
% Adam Baker 2014


%%%%%%%%%%%%%%%%%%%%%%%   P A R S E   I N P U T S   %%%%%%%%%%%%%%%%%%%%%%%

% Check SPM File Specification:
if isa(S.D,'meeg')
    D = S.D;
    S.D = D.fullfile;
else
    D = spm_eeg_load(S.D);
end

% Check Headmodel Specification:
try
    S = ft_checkopt(S,'forward_meg','char',{'Single Shell','MEG Local Spheres'});
catch
    warning('MEG Forward model specification not recognised or incorrect, assigning default: "Single Shell"')
    S = ft_setopt(S,'forward_meg','Single Shell');
end

% Check Headmodel Specification:
try
    S = ft_checkopt(S,'forward_eeg','char',{'EEG BEM'});
catch
    if isfield(S,'forward_eeg')
        error('EEG Forward model specification not recognised or incorrect')
    end
end

% Check RHINO Specification:
try
    S = ft_checkopt(S,'use_rhino',{'double','logical'});
catch
    S = ft_setopt(S,'use_rhino',0);
end

% Check Structural Specification:
S.mri = char(S.mri);
if ~isempty(S.mri)
    assert(osl_util.isfile(S.mri),'Structural MRI %s not found',S.mri);
    [~,~,ext] = fileparts(S.mri);   
    if strcmp(ext,'.gz') && S.use_rhino == 0
        % Try unzipping
        fnames = gunzip(S.mri);
        assert(length(fnames)==1,'Specified .gz file contained more than one file');
        S.mri = fnames{1};
        [~,~,ext] = fileparts(S.mri);   
    end
        
    if ~strcmp(ext,'.nii') && S.use_rhino == 0
        error('S.mri must be .nii (not .nii.gz) when using SPM coregistration')
    end
end

% Check Headshape Specification:
try
    S = ft_checkopt(S,'useheadshape',{'single','double','logical'},{0,1});
catch
    warning('Headshape specification not recognised or incorrect, assigning default: "1"')
    S = ft_setopt(S, 'useheadshape', 1);
end


if S.use_rhino
    
%%%%%   R U N   C O R E G I S T R A T I O N   U S I N G   R H I N O   %%%%%
    
    S_coreg = S;
    S_coreg.modality = {};
    if isfield(S, 'forward_meg'),
        S_coreg = rmfield(S_coreg, 'forward_meg');
        S_coreg.modality(end+1) = {'MEG'};
    end
    if isfield(S, 'forward_eeg'),
        S_coreg = rmfield(S_coreg, 'forward_eeg');
        S_coreg.modality(end+1) = {'EEG'};
    end
    S_coreg.do_plots = 0;
    D = rhino(S_coreg);
    close all
    
    
%%%%%%%%%%%%%%%%%%   R U N   F O R W A R D   M O D E L   %%%%%%%%%%%%%%%%%%
    
    D = osl_forward_model(D,S);
    D.save()
    
else % ~S.use_rhino
    
%%%%%%%   R U N   C O R E G I S T R A T I O N   U S I N G   S P M   %%%%%%%
    
    matlabbatch{1}.spm.meeg.source.headmodel.D       = {S.D};
    matlabbatch{1}.spm.meeg.source.headmodel.val     = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    
    if isempty(S.mri)
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {''};
    else
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {[S.mri ',1']};
    end;
       
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres    = 2;
    
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = S.fid.label.nasion;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = S.fid.label.lpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = S.fid.label.rpa;
    
    if(isfield(S,'fid_mnicoords'))
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = S.fid.mnicoords.nasion;
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = S.fid.mnicoords.lpa;
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = S.fid.mnicoords.rpa;
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

D = spm_eeg_load(S.D); % Reload the file

end%osl_headmodel


