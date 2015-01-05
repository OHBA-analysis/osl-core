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


% --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 


% %%%%%%%%%%%%%%%%%%   S E T   U P   T R A   M A T R I X   %%%%%%%%%%%%%%%%%%
% % This must be done for MEG Local Spheres forward model. This is because
% % ft_headmodel_localspheres.m matches the magnetometers/gradiometers to the
% % channels by looking for nonzero values in the tra matrix. Hence any
% % alterations to the tra matrix, such as Africa, or other montages must be
% % undone prior to running the forward model. As a workaround we run the
% % forward model on the pre-montaged data and copy this to the current SPM
% % MEEG object
%     
% if strcmp(S.forward_meg,'MEG Local Spheres') && isfield(D,'raw_sensors') 
%     D = spm_eeg_load(S.D);
% 
%     % check local spheres model and that raw sensor information is present
%     assert(strcmp(S.forward_meg,'MEG Local Spheres') && isfield(D,'raw_sensors'), ...
%         ['Expecting to operate on Local Spheres model, ',                      ...
%         'with raw_sensors information provided. \n']);
%     
%     % Set up a temporary SPM object to run through the matlabbatch
%     S_copy.D = D;
%     switch(lower(spm('ver')))
%         case 'spm8'
%             S_copy.newname = fullfile(D.path, ['temp_' D.fname]);
%         case {'spm12','spm12b'}
%             S_copy.outfile = fullfile(D.path, ['temp_' D.fname]);
%     end
%     Dnew = spm_eeg_copy(S_copy);
%     
%     % get original sensor object
%     sens = D.sensors('MEG');
%     
%     % ensure all labels match up and sizes are the same
%     if ~isempty(setxor(sens.label, D.raw_sensors.label)),
%         error([mfilename ':LabelMismatch'], ...
%             'Label inconsistency between D.raw_sensors and sensors(D). \n');
%     end
%     if ~isequal(size(sens.tra), size(D.raw_sensors.tra)),
%         error([mfilename ':TraSizeMismatch'], ...
%             'Size inconsistency between D.raw_sensors and sensors(D). \n');
%     end
%     
%     % set new tra matrix and match labels
%     [~,newInd] = ismember(sens.label, D.raw_sensors.label);
%     sens.tra     = D.raw_sensors.tra(newInd, :);
%     
%     Dnew = sensors(Dnew,'MEG',sens);
%     
%     Dnew.save;
%     
%     S.D = fullfile(Dnew.path,Dnew.fname);
% 
% elseif strcmpi(S.forward_meg, 'MEG Local Spheres')
%     error([mfilename ':Raw_sensorsNotSpecified'], ...
%           ['D.raw_sensors not detected. \n',      ...
%            'D.raw_sensors must be specified to allow matching of channel labels between sensor objects. \n']);
% end%if

% --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 



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



% --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 

% 
% %%%%%%%%%%   C O P Y   T O   O R I G I N A L   S P M   F I L E   %%%%%%%%%%
% 
% if exist('Dnew','var') == 1 % MEG Local Spheres case
%     Dnew  = spm_eeg_load(S.D); % reload temporary SPM object with new forward model
%     D.inv = Dnew.inv;
%     D.inv{1}.datareg.sensors.tra = D.sensors('MEG').tra;
%     D.save;
%     Dnew.delete;
% else
%     D = spm_eeg_load(S.D); % reload SPM object with new forward model
% end

% --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 8< --- 


end%osl_headmodel


