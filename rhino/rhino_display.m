function rhino_display(D,hf)
% Display of the coregistered RHINO meshes and polhemus/sensor locations.
%
% RHINO_DISPLAY(D)    - displays the coregistration for MEEG object D in a 
%                       new figure
%
% RHINO_DISPLAY(D,hf) - displays the coregistration for MEEG object D in an
%                       existing figure with handle hf
% 
% Adam Baker 2015


try 
    D.inv{1}.mesh.tess_rhino;
catch
    error('RHINO has not been run!')
end




if nargin < 2
    hf = figure;
end

if any(strcmp(get(get(hf,'children'),'type'),'axes'))
    [az,el] = view;
else
    az = 125; el = 15;
end

ha = axes('parent',hf); hold on

rotate3d(hf,'ON');

D = montage(D,'switch');

mesh       = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);

% Reduce mesh size for faster plotting
mesh_rhino = struct('faces',mesh_rhino.face,'vertices',mesh_rhino.vert);
%reduce_factor = 1e3/size(mesh_rhino.vertices,1);
%mesh_rhino = reducepatch(mesh_rhino,reduce_factor);

headshape_polhemus = D.inv{1}.datareg(1).fid_eeg.pnt;
fid_polhemus = D.inv{1}.datareg(1).fid_eeg.fid.pnt;
fid_mri = D.inv{1}.datareg(1).fid_mri.fid.pnt;

patch(mesh_rhino,'FaceColor',[238,206,179]./255,'EdgeColor','none','FaceAlpha',0.6,'Parent',ha);
patch(struct('faces',mesh.tess_ctx.face,'vertices',mesh.tess_ctx.vert),'FaceColor',[1 0.7 0.7],'EdgeColor','none','Parent',ha);
%patch(struct('faces',mesh.tess_scalp.face,'vertices',mesh.tess_scalp.vert),'EdgeColor','g','FaceColor','none','Parent',ha);

view([az,el]);

axis(ha,'image','off')
material shiny
lighting gouraud
hl = camlight('headlight');
set(hf,'WindowButtonMotionFcn',@(~,~) camlight(hl,'headlight'),'Interruptible','off','busyaction','cancel');


% Headshape locations
%--------------------------------------------------------------------------
if ~isempty(headshape_polhemus)
    %h_hsp   = plot3(headshape_polhemus(:,1),headshape_polhemus(:,2),headshape_polhemus(:,3),'dm');
    %set(h_hsp,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','r');

    err = rhino_cperror(mesh_rhino.vertices,headshape_polhemus);
    err(err>5) = 5;
    h_hsp = scatter3(headshape_polhemus(:,1),headshape_polhemus(:,2),headshape_polhemus(:,3),50,err,'filled');
    colormap(jet)
end

% Sensors (coreg.)
%--------------------------------------------------------------------------
h_sens  = plot3(D.sensors('MEG').chanpos(:,1),D.sensors('MEG').chanpos(:,2),D.sensors('MEG').chanpos(:,3),'og');
set(h_sens,'MarkerFaceColor','g','MarkerSize', 12,'MarkerEdgeColor','k');

% EEG fiducials or MEG coils (coreg.)
%--------------------------------------------------------------------------
h_fid   = plot3(fid_polhemus(:,1),fid_polhemus(:,2),fid_polhemus(:,3),'oc');
set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');

% MRI fiducials
%--------------------------------------------------------------------------
h_fidmr = plot3(fid_mri(:,1),fid_mri(:,2),fid_mri(:,3),'dm');
set(h_fidmr,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');

material dull % Improve visibility of the brain by reducing reflection glare



