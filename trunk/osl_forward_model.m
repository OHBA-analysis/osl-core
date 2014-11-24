function D = osl_forward_model(S)
%OSL_FORWARD_MODEL runs MEG or EEG forward model
%
% D = osl_forward_model(S)
% 
% Required inputs:
% S.D: SPM MEG object filename
% S.mri: structural MRI nii file name (set S.mri=[] or '' to use template
% structural)
% S.useheadshape: set to 0 or 1 to indicated if the headshape points should
% be used in the registration
%
% Optional:
% S.forward_meg: e.g. set to 'Single Shell' or 'MEG Local Spheres' (default)
%
% S.fid_label.nasion: default: 'Nasion'
% S.fid_label.lpa:    default: 'LPA';
% S.fid_label.rpa:    default: 'RPA';
%
% S.fid_mnicoords: Specify fiducual MNI coordinates with fields:
%                     .nasion - [1 x 3]
%                     .lpa    - [1 x 3]
%                     .rpa    - [1 x 3]
%                   (omit field to use SPM defaults)
%
%
% Check the output using:
% spm_eeg_inv_checkmeshes(D);
% spm_eeg_inv_checkdatareg(D);
% spm_eeg_inv_checkforward(D, val);
%
% MWW 2012

%% Allow branching into Adam's new coregistration code 
% This was the easiest way for me to wrap my code while retaining
% compatibility with the SPM8 coregistration. For OSL2/SPM12 I will adapt
% osl_headmodel to switch between SPM/my coregistration & this function
% will become obsolete
if isfield(S,'use_rhino') && S.use_rhino == 1,
    D = osl_headmodel(S);
    return
end%if

    
%%%%%%%%%%%%% BELOW USED ONLY IF NOT ROMPING WITH RHINO %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs
if ~isfield(S, 'useheadshape'), error('S.useheadshape not specified'); end

if isfield(S, 'D'),   spmfilename = S.D; else error('S.D not specified'); end
if isfield(S, 'mri'), mri = S.mri;       else error('S.mri not specified'); end

if isfield(S, 'forward_meg'), 
    forward_meg = S.forward_meg;   
else
    warning('Using Single Shell MEG head model'); 
    forward_meg = 'Single Shell'; 
end%if
if isfield(S, 'fid_label'),
    fid_label = S.fid_label;     
else
    warning('Using default Neuromag fid labels'); 
    fid_label.nasion = 'Nasion'; 
    fid_label.lpa    = 'LPA'; 
    fid_label.rpa    = 'RPA'; 
end % To see what these should be look at: D.fiducials; fidnew.fid.label

% Check Structural Specification:
try
    S.mri = char(S.mri);
    [pathstr,filestr,ext] = fileparts(S.mri);
    if isempty(ext) % force .nii suffix
        ext = 'nii';
    end
    if ~strcmpi(ext,'.nii') % catch .nii.gz - causes error in SPM
        error('Structural MRI should be a .nii file (not .nii.gz)')
    else
        mri = fullfile(pathstr,[filestr '.nii']); % force .nii suffix
    end
catch
    mri = ft_getopt(S,'mri','');
    warning('Structural MRI specification not recognised or incorrect, using standard structural');
end


%% Load data
D = spm_eeg_load(spmfilename);

spmfilename_orig = spmfilename; % save in case changed in next section

%% Update tra matrix for MEG local spheres
if isfield(D,'raw_sensors') && strcmpi(forward_meg, 'MEG Local Spheres');
    %%%%%%%%%%%%%%%
    % replace the tra matrix with that held from before the montage
    % This is to ensure consistency after montages are applied - e.g. from
    % normalisation, or AFRICA. Talk to Mark Woolrich
    % (mark.woolrich@ohba.ox.ac.uk) for more information
    %
    % We will need to ensure that the rows of the tra matrix match the
    % labels held in grad.label
    %
    % Updated by GC 04-12-2013
    
    % check that the old syntax is not being used
    if isfield(D, 'raw_tra_matrix'),
        if isequal(D.raw_tra_matrix, D.raw_sensors.tra),
            warning([mfilename ':DeprecatedField'],               ...
                    ['Deprecated field D.raw_tra_matrix found. ', ...
                     'This will be ignored. \n']);
        else
            warning([mfilename ':DeprecatedFieldDiscrepancy'],               ...
                    ['Deprecated field D.raw_tra_matrix found. \n',          ...
                     'It will be ignored, even though its contents differ ', ...
                     'from the currently saved tra matrix. \n']);
        end%if
    else % does not exist
        % we use the tra matrix from within raw_sensors
    end%if
           
    % set up a temporary D object to run through the matlabbatch
    S_copy.D       = D;
    S_copy.newname = fullfile(D.path, ['temp_' D.fname]);
    Dnew           = spm_eeg_copy(S_copy);
    
    % get original sensor object
    grad = D.sensors('MEG');
    
    % ensure all labels match up and sizes are the same
    if ~isempty(setxor(grad.label, D.raw_sensors.label)),
        error([mfilename ':LabelMismatch'], ...
              'Label inconsistency between D.raw_sensors and sensors(D). \n');
    end%if
    if ~isequal(size(grad.tra), size(D.raw_sensors.tra)),
        error([mfilename ':TraSizeMismatch'], ...
              'Size inconsistency between D.raw_sensors and sensors(D). \n');
    end%if
    
    % set new tra matrix and permute rows to match labels already in grad
    [~, newInd] = ismember(grad.label, D.raw_sensors.label);
    new_tra     = D.raw_sensors.tra(newInd, :);
    
    grad.tra = new_tra;
    Dnew     = sensors(Dnew,'MEG',grad);  
    
    Dnew.save;
    
    % update filename passed to matlab batch
    spmfilename  = fullfile(Dnew.path, Dnew.fname);
    copy_fm_flag = 1; % flag up that we have maniuplated tra matrix

elseif strcmpi(forward_meg, 'MEG Local Spheres')
    error([mfilename ':Raw_sensorsNotSpecified'],                   ...
          ['D.raw_sensors not detected. \n',                        ...
           'D.raw_sensors must be specified to allow matching of ', ...
           'channel labels between sensor objects. \n']);
else
    % not MEG local spheres - carry on as normal
    copy_fm_flag = 0;
end%if MEG local spheres

%% Set up matlab batch
matlabbatch{1}.spm.meeg.source.headmodel.D       = {spmfilename};
matlabbatch{1}.spm.meeg.source.headmodel.val     = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';

if ~strcmp(mri,''),
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {[mri ',1']};
else
    % should not be here - set above
    warning([mfilename ':MRInotSet'], ...
          'Ideally, need an mri file set to complete the head model. \n');
end%if mri specified

matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;

matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = fid_label.nasion;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = fid_label.lpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = fid_label.rpa;
 
if(isfield(S,'fid_mnicoords'))
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = S.fid_mnicoords.nasion;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = S.fid_mnicoords.lpa;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = S.fid_mnicoords.rpa;
else
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
end%if fid_mnicoords

matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = S.useheadshape;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = forward_meg;

%% Run forward model
spm_jobman('run', matlabbatch);

%% Tidy output
% if MEG local spheres, hold the sensor object used in forward model in inv
if copy_fm_flag
    Dnew     = spm_eeg_load(spmfilename);
    D.inv{1} = Dnew.inv{1};
    D.save;
    Dnew.delete;
else
    D = spm_eeg_load(spmfilename);
end

% display recognised data type
chant = D.chantype;

if sum(ismember(unique(chant),'MEGGRAD')), 
    disp('Using CTF data.');
    
elseif sum(ismember(unique(chant),'MEGPLANAR')) && ~sum(ismember(unique(chant),'MEGPCACOMP'))
    disp('Using Elekta Neuromag 306 data.');       
    
elseif sum(ismember(unique(chant),'EEG')), % added by DM
    disp('Using EEG data.');
    
else
    warning([mfilename ': Unrecognised data source. \n']);
end%if
end%osl_forward_model
