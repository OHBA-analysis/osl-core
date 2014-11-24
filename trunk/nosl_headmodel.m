function D = nosl_headmodel(S)
% NOSL_HEADMODEL runs MEG coregistration and forward model in SPM8 or 
% SPM12. These two tasks are separately performed by nosl_datareg.m and 
% nosl_forward_model.m which are wrapped together to ensure functionality 
% of local spheres forward model when montaging has been applied
%
% D = nosl_headmodel(S)
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
% S.neuromag_planar_baseline_correction: determines whether this correction
%                                        is performed. 
%
% Adam Baker 2014




%%%%%%%%%%%%%%%%%%%%%%%   P A R S E   I N P U T S   %%%%%%%%%%%%%%%%%%%%%%%

% Check SPM File Specification:
try
    S.D = char(S.D);
    [pathstr,filestr] = fileparts(S.D);
    S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
    D = spm_eeg_load(S.D);
catch ME
    error(ME.identifier, ...
          'SPM file specification not recognised or incorrect\n');
end%try

% Check Headmodel Specification:
try
    S = ft_checkopt(S,'forward_meg','char',{'Single Shell','MEG Local Spheres'});
catch ME
    warning(ME.identifier,                                                ...
            ['Forward model specification not recognised or incorrect, ', ...
             'assigning default: "Single Shell"'])
    S = ft_setopt(S,'forward_meg','Single Shell');
end%try


%%%%%%%%%%%%%%%%%%   S E T   U P   T R A   M A T R I X   %%%%%%%%%%%%%%%%%%
% This must be done for MEG Local Spheres forward model. This is because
% ft_headmodel_localspheres.m matches the magnetometers/gradiometers to the
% channels by looking for nonzero values in the tra matrix. Hence any
% alterations to the tra matrix, such as Africa, or other montages must be
% undone prior to running the forward model. As a workaround we run the
% forward model on the pre-montaged data and copy this to the current SPM
% MEEG object
    
if strcmp(S.forward_meg,'MEG Local Spheres') && isfield(D,'raw_sensors') 
    S = set_up_tra(S);
    
elseif strcmpi(S.forward_meg, 'MEG Local Spheres')
    error([mfilename ':Raw_sensorsNotSpecified'], ...
          ['D.raw_sensors not detected. \n',      ...
           'D.raw_sensors must be specified to allow matching of channel labels between sensor objects. \n']);
end%if


%%%%%   R U N   C O R E G I S T R A T I O N   U S I N G   R H I N O   %%%%%

% input parsing done within Rhino. 
% strip out un-necessary fields
S_coreg = S;
if isfield(S_coreg, 'neuromag_planar_baseline_correction'),
    S_coreg = rmfield(S_coreg, 'neuromag_planar_baseline_correction');
end%if
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
nosl_forward_model(S_forward);


%%%%%%%%%%   C O P Y   T O   O R I G I N A L   S P M   F I L E   %%%%%%%%%%

if exist('Dnew','var') == 1 % MEG Local Spheres case
    Dnew  = spm_eeg_load(S.D); % reload temporary SPM object with new forward model
    D.inv = Dnew.inv;
    D.save;
    Dnew.delete;
else
    D = spm_eeg_load(S.D); % reload SPM object with new forward model
end

    
%%%%%%%   N E U R O M A G   B A S E L I N E   C O R R E C T I O N   %%%%%%%
if  sum(ismember(unique(D.chantype),'MEGPLANAR')) && ...
   ~sum(ismember(unique(D.chantype),'MEGPCACOMP')),

    D = apply_neuromag_baseline_correction(S);
end%if
end%nosl_headmodel

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = set_up_tra(S)
%SET_UP_TRA sets up tra matrix for local spheres model
% This must be done for MEG Local Spheres forward model. This is because
% ft_headmodel_localspheres.m matches the magnetometers/gradiometers to the
% channels by looking for nonzero values in the tra matrix. Hence any
% alterations to the tra matrix, such as Africa, or other montages must be
% undone prior to running the forward model. As a workaround we run the
% forward model on the pre-montaged data and copy this to the current SPM
% MEEG object

D = spm_eeg_load(S.D);

% check local spheres model and that raw sensor information is present
assert(strcmp(S.forward_meg,'MEG Local Spheres') && isfield(D,'raw_sensors'), ...
       ['Expecting to operate on Local Spheres model, ',                      ...
        'with raw_sensors information provided. \n']);

% Set up a temporary SPM object to run through the matlabbatch
S_copy.D = D;
switch(lower(spm('ver')))
    case 'spm8'
        S_copy.newname = fullfile(D.path, ['temp_' D.fname]);
    case {'spm12','spm12b'}
        S_copy.outfile = fullfile(D.path, ['temp_' D.fname]);
end
Dnew = spm_eeg_copy(S_copy);

% get original sensor object
sens = D.sensors('MEG');

% ensure all labels match up and sizes are the same
if ~isempty(setxor(sens.label, D.raw_sensors.label)),
    error([mfilename ':LabelMismatch'], ...
        'Label inconsistency between D.raw_sensors and sensors(D). \n');
end
if ~isequal(size(sens.tra), size(D.raw_sensors.tra)),
    error([mfilename ':TraSizeMismatch'], ...
        'Size inconsistency between D.raw_sensors and sensors(D). \n');
end

% set new tra matrix and match labels
[~,newInd] = ismember(sens.label, D.raw_sensors.label);
sens.tra     = D.raw_sensors.tra(newInd, :);

Dnew = sensors(Dnew,'MEG',sens);

Dnew.save;

% update filename passed to matlab batch
S.D = fullfile(Dnew.path,Dnew.fname);
end%set_up_tra
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = apply_neuromag_baseline_correction(S)
    
    D = spm_eeg_load(S.D);
    
    % check neuromag
    assert( sum(ismember(unique(D.chantype),'MEGPLANAR')) && ...
           ~sum(ismember(unique(D.chantype),'MEGPCACOMP')),  ...
           'Expecting Elekta Neuromag data. \n');
       
    disp('Using Elekta Neuromag 306 data.');
    
    % check to see if baseline correction has been done before or not
    grad = D.sensors('MEG');
    if ~isfield(grad,'tra_baseline_corrected')
        grad.tra_baseline_corrected = 0;
    end%if
    
    if grad.tra_baseline_corrected
        % if done before then do not redo now:
        if isfield(S,'neuromag_planar_baseline_correction') && ~strcmp(S.neuromag_planar_baseline_correction,'none')
            warning('No gradiometer baseline correction being carried out, as it has already been done.');
        end%if
        
    else
        if  isfield(S,'neuromag_planar_baseline_correction') && ~strcmp(S.neuromag_planar_baseline_correction,'none'),
            % do baseline correction
            D = osl_neuromag_grad_baseline_correction(S.D, S.neuromag_planar_baseline_correction);
            
        else
            % baseline correction not been done, and is not being done now
            if grad.tra_baseline_corrected == 0,
                warning(['You have not gradiometer baseline corrected - you may need to make sure this is done. Note '...
                         'that this should be done before any montaging (e.g. before AFRICA), and can be done by turning '...
                         'on the opt.coreg.neuromag_planar_baseline_correction=1 option if using OPT.)']);
            end%if            
        end%if
    end%if
end%apply_neuromag_baseline_correction
% [EOF]
