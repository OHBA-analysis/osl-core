function D = rhino(S)
% RHINO - Registration using Headshapes Including Nose in OSL
% performs co-registration of two sets of fiducials according 
% to sets of corresponding points and headshape points. Uses FLIRT & 
% BETSURF for scalp extraction and multistart ICP for headshape matching.
% D = rhino(S)
%
% REQUIRED INPUTS:
%
% S.D                - SPM MEG object filename or MEEG object
%
% S.mri              - structural MRI .nii/.nii.gz file name OR [] or '' to use template
%
%
% OPTIONAL INPUTS 
%
% S.do_plots         - produce diagnostic plots 
%                        (default = 0)
%
% S.multistart       - number of random initialisations of headshape fitting
%                        (default 10)
%
% S.useheadshape     - used headshape points for the coregistration (0 or 1)
%                        (default = 1)
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
% S.modality          - Modalities to use {'MEG','EEG'} (default, 'MEG')
%
%
% S.fid_headcastcoords  - Specify headcast native coordinates with fields:
%                           .LB - [1 x 3]
%                           .LF - [1 x 3]
%                           .RB - [1 x 3]
%                           .RF - [1 x 3]
%
%
%                            ,-.             __
%                          ,'   `---.___.---'  `.
%                        ,'   ,-                 `-._
%                      ,'    /                       \
%                   ,\/     /                        \\
%               )`._)>)     |                         \\
%               `>,'    _   \                  /       ||
%                 )      \   |   |            |        |\\
%        .   ,   /        \  |    `.          |        | ))
%        \`. \`-'          )-|      `.        |        /((
%         \ `-`   .`     _/  \ _     )`-.___.--\      /  `'
%          `._         ,'     `j`.__/           `.    \
%            / ,    ,'         \   /`             \   /
%            \__   /           _) (               _) (
%              `--'           /____\             /____\
%
% Adam Baker 2014


% Todo
% - potentially split code to separate scalp extraction & coregistration


%%%%%%%%%%%%%%%%%%%%%%%   P A R S E   I N P U T S   %%%%%%%%%%%%%%%%%%%%%%%

disp('*** RUNNING RHINO COREGISTRATION ***')

% Check SPM File Specification:
if isa(S.D,'meeg')
    D = S.D;
    S.D = D.fullfile;
else
    D = spm_eeg_load(S.D);
end

% Check Structural Specification:
S.mri = char(S.mri);
if isempty(S.mri)
    tempMesh = spm_eeg_inv_mesh;
    sMRI     = tempMesh.sMRI;
    warning('Using SPM template MRI');
else
    [~,~,ext] = fileparts(S.mri);
    if strcmp(ext,'.nii')
        sMRI = S.mri;
    else
        fnames = gunzip(S.mri);
        assert(length(fnames)==1,'Specified .gz file contained more than one file');
        S.mri = fnames{1}; % Swap in the nii file
        sMRI = S.mri;
    end
end
[~,~,ext] = fileparts(sMRI);   
assert(strcmp(ext,'.nii'),'Unexpected extension - perhaps the .gz file did not contain a .nii file?'); % We can assume after this point we are working with a .nii file
assert(osl_util.isfile(sMRI),'MRI file %s not found',sMRI)

% Check Headshape Specification:
try
    S = ft_checkopt(S,'useheadshape',{'single','double','logical'},{0,1});
catch
    warning('Headshape specification not recognised or incorrect, assigning default: "1"')
    S = ft_setopt(S, 'useheadshape', 1);
end%try

% Check headshape points exist
if S.useheadshape == 1 && isempty(D.fiducials.pnt)
    error('SPM file doesn''t contain any headshape points!')
end

% Check Headshape Specification:
if ~isfield(S,'modality') || isempty(S.modality)  
    S.modality = {'MEG'};
end

% Check Fiducial Specification:
if ~isfield(S,'fid') || isempty(S.fid)  
    S.fid = struct;
    manual_fid_selection = 1;
else
    manual_fid_selection = 0;
end

try
    S = ft_checkopt(S,'fid','struct');
catch
    error('Fiducial specification should be a structure. \n');
end

if ~isfield(S.fid,'coordsys') || isempty(S.fid.coordsys)
    S.fid.coordsys = 'MNI';
end

if ~isfield(S.fid,'coords') || isempty(S.fid.coords)
    if ~strcmpi(S.fid.coordsys,'mni')
        warning('No coordinates specified, using default MNI coordinates')
        S.fid.coordsys = 'mni';
    end
    S.fid.coords.nasion = [  1   85  -41];
    S.fid.coords.rpa    = [ 83  -20  -65];
    S.fid.coords.lpa    = [-83  -20  -65];
else
    try
        S.fid = ft_checkopt(S.fid,'coords','struct');
    catch
        error('Fiducial coordinate specification should be a structure. \n');
    end
end

if ~strcmpi(S.fid.coordsys,'headcast')
    try
        S.fid = ft_checkopt(S.fid,'label','struct');
        assert(isfield(S.fid.label, 'nasion') &&        ...
            isfield(S.fid.label, 'lpa')    &&        ...
            isfield(S.fid.label, 'rpa'),             ...
            [mfilename ':fid.labelIncorrectFields'], ...
            'Incorrect fields in S.fid.label\n');
    catch
        warning('Fiducial label specification not recognised or incorrect, assigning Elekta defaults\n')
        S.fid = ft_setopt(S.fid,'label',struct('nasion','Nasion','lpa','LPA','rpa','RPA'));
    end
    
end

if numel(fieldnames(S.fid.label)) == 3
    fid_labels = {S.fid.label.nasion, S.fid.label.rpa, S.fid.label.lpa};
else
    fid_labels = fieldnames(S.fid.label)';
end

switch lower(S.fid.coordsys)
    case 'native'
       fid_native  = [S.fid.coords.nasion; ...
                      S.fid.coords.rpa;    ...
                      S.fid.coords.lpa];
    case 'mni'
       fid_MNI     = [S.fid.coords.nasion; ...
                      S.fid.coords.rpa;    ...
                      S.fid.coords.lpa];
    case 'headcast'
       fid_native  = [S.fid.coords.LB;  ...
                      S.fid.coords.LF;  ...
                      S.fid.coords.RB;  ...
                      S.fid.coords.RF];
    otherwise
        error('S.fid.coordsys must be either ''Native'' or ''MNI''')
end


% Check Plotting Specification:
S.do_plots   = ft_getopt(S,'do_plots',   0);
S.multistart = ft_getopt(S,'multistart', 10);


%%%%%%%%%%%%%%%%   G E T   S C A L P   U S I N G   F S L   %%%%%%%%%%%%%%%%

% Copy sMRI to new file for modification
[struct_path, original_name, ext] = fileparts(sMRI); 
struct_name = sprintf('rhino_%s',original_name);
copyfile(fullfile(struct_path,[original_name ext]),fullfile(struct_path,[struct_name ext]));
sMRI = fullfile(struct_path,[struct_name ext]); % Point to the copied file

% File names (generated by FSL calls)
mesh_file  = fullfile(struct_path,[struct_name '_mesh.vtk']);
trans_file = fullfile(struct_path,[struct_name '_trans.txt']);
bet_file   = fullfile(struct_path,[struct_name '_outskin_mesh.nii.gz']);
scalp_file = fullfile(struct_path,[struct_name '_scalp.nii.gz']);
std_brain  = [getenv('FSLDIR') '/data/standard/MNI152_T1_1mm.nii.gz'];

% determine orientation of standard brain
cmd = sprintf('fslorient -getorient %s', std_brain);
std_orient = deblank(runcmd(cmd));
if ~ismember(std_orient,{'RADIOLOGICAL','NEUROLOGICAL'})
    error('Cannot determine orientation of standard brain, please check output of:\n%s',cmd)
end

% determine orientation of subject brain
cmd = sprintf('fslorient -getorient %s', sMRI);
smri_orient = deblank(runcmd(cmd));
if ~ismember(smri_orient,{'RADIOLOGICAL','NEUROLOGICAL'})
    error('Cannot determine orientation of subject brain, please check output of:\n%s',cmd);
end

% if orientations don't match, reorient the subject brain
if ~strcmp(std_orient, smri_orient)
    disp('reorienting subject brain to match standard brain');
    switch std_orient
        case 'RADIOLOGICAL'
            orient_flag = '-forceradiological';
        case 'NEUROLOGICAL'
            orient_flag = '-forceneurological';
        otherwise
            error('cannot determine orientation of subject brain');
    end
    runcmd('fslorient %s %s', orient_flag, sMRI);
end

% CHECK IF SCALP EXTRACTION ALREADY DONE
if ~osl_util.isfile(scalp_file)
    
    % RUN FLIRT - get the transformation to take the input image into standard MNI space
    if ~osl_util.isfile(trans_file)
        disp('Running FLIRT...')
        flirtCommand = runcmd('flirt -in %s -ref %s -omat %s -o %s',sMRI, std_brain, trans_file,fullfile(struct_path,[struct_name,'_MNI']));
    end
    
    % RUN BET - get the surface mesh in the scanner space
    if ~osl_util.isfile(mesh_file)
        disp('Running BET...')
        runcmd('bet2 %s %s -n --mesh ',sMRI,fullfile(struct_path,struct_name)); % This generates the mesh file
    end
    
    % RUN BETSURF - get the head surfaces in scanner space, using FLIRTs trans_file as part of the extraction (see betsurf documentation)
    if ~osl_util.isfile(bet_file)
        disp('Running BETSURF...')
        runcmd('betsurf --t1only -o %s %s %s %s',sMRI,mesh_file,trans_file,fullfile(struct_path,struct_name))
    end
    
    % CLEAN UP
    % Use system instead of dos because an error might occur here
    delete([fullfile(struct_path,struct_name) '*.off']);
    delete([fullfile(struct_path,struct_name) '_outskull*.*']);
    delete([fullfile(struct_path,struct_name) '_inskull*.*']);
    delete([mesh_file]);
    

    %%%%%%%%%%%%%%%   R E F I N E   S C A L P   O U T L I N E   %%%%%%%%%%%%%%%

    disp('Running scalp extraction')

    % READ IN VOLUME & OUTLINE
    
    % check it has been created properly
    assert( osl_util.isfile(bet_file),    ...
           [mfilename ':BET_FILEDoesNotExist'], ...
           ['bet_file does not exist. '         ...
            'Maybe it failed to create earlier in Rhino?\n']);
       
    scalp            = read_avw(bet_file);
    [vol, ~, scales] = read_avw(sMRI);
    
    vol = vol(:,:,:,1); % In case of 4D volumes
    vol = vol ./ max(vol(:));
    
    % PLOT VOLUME AND BETSURF OUTLINE:
    if S.do_plots
        hf = figure; ha = axes('parent',hf);
        for i = 1:size(vol,3)
            im = repmat(vol(:,:,i),[1,1,3]);
            im_outline = zeros(size(im));
            im_outline(:,:,1) = scalp(:,:,i);
            im(repmat(im_outline(:,:,1),[1,1,3])==1) = 0;
            im = im + im_outline;
            image(im,'parent',ha)
            axis(ha,'image','off')
            drawnow
        end
        close(hf)
    end
    
    % CREATE MASK BY FILLING OUTLINE
    mask = ones(size(scalp)+2); % add border to account for gaps at edges
    mask(2:end-1,2:end-1,2:end-1) = scalp;
    mask([1:3,end-2:end],[1:3,end-2:end],[1:3,end-2:end]) = 0;
    mask = imfill(mask,'holes');
    mask = mask(2:end-1,2:end-1,2:end-1);
    
    % RECLASSIFY BRIGHT VOXELS OUTSIDE OF MASK (E.G. NOSE)
    vol_inside = vol;  vol_inside( mask==0) = nan;
    vol_outside = vol; vol_outside(mask==1) = nan;
    
    mix = gmdistribution(cat(1,nanmean(vol_outside(:)),nanmean(vol_inside(:))), ...
                         cat(3,nanvar(vol_outside(:)),nanvar(vol_inside(:))), ...
                         cat(1,sum(~isnan(vol_outside(:))),sum(~isnan(vol_inside(:)))));
    
    vol_outside(isnan(vol_outside)) = 0;
    idx = cluster(mix,vol_outside(:));
    mask_class = zeros(size(vol));
    mask_class(1:size(idx)) = idx;
    
    % MAKE NEW MASK FROM ORIGINAL AND RECLASSIFIED VOXELS
    mask = double(mask_class==2 | mask==1);
    
    % CLEAN UP MASK
    mask = imfill(mask,'holes');
    for i = 1:size(mask,1)
        mask(i,:,:) = bwmorph(squeeze(mask(i,:,:)),'majority');
    end
    for j = 1:size(mask,2)
        mask(:,j,:) = bwmorph(squeeze(mask(:,j,:)),'majority');
    end
    for k = 1:size(mask,3)
        mask(:,:,k) = bwmorph(squeeze(mask(:,:,k)),'majority');
    end
    mask = imfill(mask,'holes');
    
    % FILL HOLES PER SLICE (USEFUL FOR FILLING NOSTRILS)
    for i = 1:size(mask,1)
        mask(i,:,:) = imfill(squeeze(mask(i,:,:)),'holes');
    end
    for j = 1:size(mask,2)
        mask(:,j,:) = imfill(squeeze(mask(:,j,:)),'holes');
    end
    for k = 1:size(mask,3)
        mask(:,:,k) = imfill(squeeze(mask(:,:,k)),'holes');
    end
    
    % EXTRACT OUTLINE
    outline = zeros(size(mask));
    for i = 1:size(mask,1)
        outline(i,:,:) = bwmorph(squeeze(mask(i,:,:)),'remove');
    end
    for j = 1:size(mask,2)
        outline(:,j,:) = bwmorph(squeeze(mask(:,j,:)),'remove');
    end
    for k = 1:size(mask,3)
        outline(:,:,k) = bwmorph(squeeze(mask(:,:,k)),'remove');
    end
    
    % SET OUTLINE AS 1, INSIDE AS 0 AND BACKGROUND AS 2:
    mask(mask==0)    = 2;
    mask(mask==1)    = 0;
    mask(outline==1) = 1;
    outline = mask;
    
    % SAVE AS NIFTI
    save_avw(outline,scalp_file,'d',scales);
    runcmd(['fslcpgeom  ' bet_file ' ' scalp_file], 1);
    
    % PLOT VOLUME AND NEW OUTLINE:
    if S.do_plots
        hf = figure; ha = axes('parent',hf);
        for i = 1:size(vol,3)
            im = repmat(vol(:,:,i),[1,1,3]);
            im_outline = zeros(size(im));
            im_outline(:,:,1) = outline(:,:,i);
            im(repmat(im_outline(:,:,1),[1,1,3])==1) = 0;
            im = im + im_outline;
            image(im,'parent',ha)
            axis(ha,'image','off')
            drawnow
        end
        close(hf)
    end
        
else
    disp(['Loading existing scalp extraction from ' scalp_file])
    outline = read_avw(scalp_file);  
end

if isfield(D,'inv')
    D = rmfield(D,'inv');
end

% COMPUTE SPM MESHES
disp('Computing SPM meshes')
D.inv{1}.mesh = spm_eeg_inv_mesh(sMRI,2);

% CREATE SURFACE MESH ON WHICH TO OVERLAY ICP FIT
mask = outline; mask(mask~=2) = 1; mask(mask==2) = 0;
%mask = permute(mask,[2 1 3]); % no idea why I need to permute this...
[x,y,z] = meshgrid(1:size(mask,1), 1:size(mask,2), 1:size(mask,3));
x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);

mask = smooth3(mask,'gaussian',[9,9,9]);
mesh_rhino = isosurface(x,y,z,mask);
mesh_rhino = reducepatch(mesh_rhino,0.1);


%%%%%%%%%%%%%%%   P E R F O R M   R E G I S T R A T I O N   %%%%%%%%%%%%%%%

disp('Coregistering')

% MAP INTO SCANNER COORDINATE SYSTEM
% Have to be really careful using qform & sform. Sometimes they are not
% valid, or the wrong code is used meaning a different tranformation is
% applied to the one used in FLIRT. If the output from RHINO looks silly 
% then the qform/sform or qformcode/sformcode are likely causes. If in
% doubt then check that the headshape points in all coordinate systems
% correspond to those in FSLVIEW.
qformcode = str2double(runcmd(['fslorient -getqformcode ' scalp_file], 1));
sformcode = str2double(runcmd(['fslorient -getsformcode ' scalp_file], 1));
% qform and/or sform need to be valid to keep orientation consistent:
if qformcode == 1 % code == 1 usually means voxels to anatomical units
    qform_native = runcmd(['fslorient -getqform ' sMRI]);
elseif sformcode == 1
    warning('qform code is not valid, using sform \n');
    qform_native = runcmd(['fslorient -getsform ' sMRI]);
elseif sformcode == 2
    warning('Neither qform or sform codes == 1, using sform: please check output carefully');
    qform_native = runcmd(['fslorient -getsform ' sMRI]);
elseif qformcode == 2
    warning('Neither qform or sform codes == 1, using qform: please check output carefully');
    qform_native = runcmd(['fslorient -getqform ' sMRI]);
else
    warning(['Neither sform or qform codes are valid... Trying anyway: ', ...
             'Please check the registration carefully for errors.\n']);
    qform_native = runcmd(['fslorient -getqform ' sMRI]);
end
qform_native = str2num(qform_native);                                      
qform_native = reshape(qform_native,4,4)';

[x,y,z] = ind2sub(size(outline),find(outline == 1));
headshape = [x y z];
headshape_native = spm_eeg_inv_transform_points(qform_native,headshape);
mesh_rhino.vertices = spm_eeg_inv_transform_points(qform_native,mesh_rhino.vertices);


% GET TRANSFORM FROM MNI TO NATIVE
qform_mni = runcmd(['fslorient -getqform ' std_brain]);
qform_mni = str2num(qform_mni);                                            
qform_mni = reshape(qform_mni,4,4)';
toMNI = load(trans_file);
% trans seems to map from native (coords) to MNI (slices) so the correct
% toMNI transformation is qform_mni*trans_file
%toMNI = qform_mni*toMNI*inv(qform_native);

% GET POLHEMUS FIDUCIALS
fid_polhemus       = D.fiducials.fid.pnt;
polhemus_labels    = D.fiducials.fid.label;
assert(~isempty(fid_polhemus)       &&  ...
       ~isempty(polhemus_labels),       ...
       [mfilename ':NoPolhemusPoints'], ...
       'Polhemus points not found in D.fiducials. \n');

% TRANSFORM MNI FIDUCIALS TO NATIVE SPACE IF DEFINED
if exist('fid_MNI','var')
    if ~isempty(fid_MNI)
        %fid_native = [fid_MNI, ones(3,1)] * inv(toMNI)';
        fid_native = [fid_MNI, ones(3,1)] * inv(D.inv{1}.mesh.Affine)';
        fid_native = fid_native(:,1:3);
    end
end

% MATCH FIDUCIALS
[~,fid_order] = ismember(arrayfun(@(x) {lower(polhemus_labels{x}(1:min(3,length(polhemus_labels{x}))))}, 1:length(fid_labels)),...
                         arrayfun(@(x) {lower(fid_labels{x}(1:min(3,length(fid_labels{x}))))}, 1:length(fid_labels)));

fid_labels = fid_labels(fid_order);
fid_native = fid_native(fid_order, :);

fid_polhemus = fid_polhemus(1:length(fid_order),:);

% REFINE FIDUCIALS MANUALLY
if manual_fid_selection
    fid_native = rhino_select(mesh_rhino,fid_native,fid_labels);
end

% RIGID BODY TRANSFORMATION USING FIDUCIALS
M_rigid         = spm_eeg_inv_rigidreg(fid_polhemus',   fid_native');
headshape_coreg = spm_eeg_inv_transform_points(M_rigid, headshape_native);
fid_coreg       = spm_eeg_inv_transform_points(M_rigid, fid_native);

% APPLY TRANSFORMATION TO SCALP MESH
mesh_rhino.vertices = spm_eeg_inv_transform_points(M_rigid,mesh_rhino.vertices);

% SET UP PLOT
if S.do_plots
    hf       = figure; 
    ha       = axes('parent',hf);
    patch(struct('faces',    mesh_rhino.faces,     ...
                 'vertices', mesh_rhino.vertices), ...
          'FaceColor', [238,206,179]./255,          ...
          'EdgeColor', 'none',                      ...
          'FaceAlpha', 0.5,                         ...
          'Parent', ha);
    view(90,0); 
    axis(ha,'image','off')
    material shiny
    lighting gouraud
    camlight(0,0);
    rotate3d(hf, 'on');
else
    hf = [];
    ha = [];
end%if
% neat cleanup, even if stopped by error
close_hf = onCleanup(@() close(hf));

% ICP RIGID BODY TRANSFORMATION USING HEADSHAPE POINTS AND SCALP OUTLINE
if S.useheadshape
    % get headshape points
    headshape_polhemus = [D.fiducials.pnt; D.fiducials.fid.pnt];
    assert(~isempty(headshape_polhemus),      ...
           [mfilename ':NoHeadhshapePoints'], ...
           'Headshape points not found in D.fiducials.pnt. \n');
    % Reset random number seed (to ensure consistency over multiple sessions)
    rng(1,'twister'); 
    % Run ICP with multiple initialisations
    M_icp_inv = rhino_icp(headshape_coreg', headshape_polhemus', S.multistart, ha);
    M_icp     = inv(M_icp_inv);
    
%     headshape_icp = spm_eeg_inv_transform_points(M_icp, headshape_coreg); 
%     Mgd = rhino_gradientdescent(headshape_icp,headshape_polhemus);
    
    % headshape_coreg = spm_eeg_inv_transform_points(M_icp,headshape_coreg);
    fid_coreg       = spm_eeg_inv_transform_points(M_icp, fid_coreg);
else
    % just keep results from rigid body
    M_icp = eye(4);
end%if

% APPLY TO SURFACE MESH VERTICES
mesh_rhino.vertices = spm_eeg_inv_transform_points(M_icp,mesh_rhino.vertices);

if S.do_plots
    clear close_hf
end

%%%%%%%%%%%%%%%   P O P U L A T E   M E E G   F I E L D S  %%%%%%%%%%%%%%%%

% MESHES IN MNI SPACE (so far haven't used scalp.mat - may be necessary)
Mgifti = M_icp * M_rigid; % Assumes SPM MNI-native transform is okay
%Mgifti = M_icp * M_rigid * inv(toMNI) * D.inv{1}.mesh.Affine; % Assumes SPM MNI-native transform is wrong

% NOTE TO MYSELF - AT THE MOMENT I AM ASSUMING THAT THE SPM MNI-NATIVE
% TRANSFORM (D.INV{1}.MESH.AFFINE) IS THE SAME AS TOMNI. IT MIGHT BE MORE
% SENSIBLE TO TRANSFROM ALL THE MESHES INTO RHINO NATIVE SPACE AND AMEND
% THE AFFINE TRANSFORMATION MATRIX AS NECESSARY.

% APPEND RHINO MESH TO SPM OBJECT
D.inv{1}.mesh.tess_rhino.face = mesh_rhino.faces;
D.inv{1}.mesh.tess_rhino.vert = mesh_rhino.vertices;

%mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);
%mesh       = spm_eeg_inv_transform_mesh(Mgifti, D.inv{1}.mesh);

scalp          = gifti(D.inv{1}.mesh.tess_scalp);
scalp.vertices = spm_eeg_inv_transform_points(Mgifti,scalp.vertices);
scalp          = struct('faces',scalp.faces,'vertices',scalp.vertices);
% ctx            = gifti(D.inv{1}.mesh.tess_ctx);
% ctx.vertices   = spm_eeg_inv_transform_points(Mgifti,ctx.vertices);
% ctx            = struct('faces',ctx.faces,'vertices',ctx.vertices);

D.inv{1}.comment = 'rhino';
D.inv{1}.date = char({date; datestr(now,'HH:MM')});

% WRITE DATAREG TO SPM OBJECT
D.inv{1}.datareg = [];
for modality = S.modality
    D.inv{1}.datareg(end+1).modality        = char(modality);
    D.inv{1}.datareg(end).sensors           = D.sensors(char(modality));
    D.inv{1}.datareg(end).fid_eeg           = D.fiducials;
    D.inv{1}.datareg(end).fid_mri.fid.label = fid_labels;
    D.inv{1}.datareg(end).fid_mri.fid.pnt   = fid_coreg;
    D.inv{1}.datareg(end).fid_mri.pnt       = scalp.vertices;
    D.inv{1}.datareg(end).fromMNI           = M_icp * M_rigid * inv(D.inv{1}.mesh.Affine);
    D.inv{1}.datareg(end).toMNI             = inv(D.inv{1}.datareg(end).fromMNI);
end

D.save;

%%%%%%%%%%%%%%%%%%%%%%%   P L O T   R E S U L T S  %%%%%%%%%%%%%%%%%%%%%%%%
if S.do_plots
    rhino_display(D);
end
disp('***  RHINO COREGISTRATION COMPLETE ***')
end%rhino
%#ok<*ST2NM> % suppress str2num warnings
%#ok<*MINV>  % suppress inv warnings as only small affine matrices are being inverted
