function rhino_display(D)

hf = figure('WindowButtonMotionFcn',@movelight,'Interruptible','off','busyaction','cancel');
ha = axes('parent',hf); hold on
rotate3d(hf);

D = montage(D,'switch');

mesh_rhino = D.inv{1}.mesh.tess_rhino;
mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg.fromMNI*D.inv{1}.mesh.Affine,mesh_rhino);

mesh = spm_eeg_inv_transform_mesh(D.inv{1}.datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
mesh_ctx = mesh.tess_ctx;

headshape_polhemus = D.inv{1}.datareg.fid_eeg.pnt;
fid_polhemus = D.inv{1}.datareg.fid_eeg.fid.pnt;
fid_mri = D.inv{1}.datareg.fid_mri.fid.pnt;

patch(struct('faces',mesh_rhino.face,'vertices',mesh_rhino.vert),'FaceColor',[238,206,179]./255,'EdgeColor','none','FaceAlpha',0.8,'Parent',ha);

patch(struct('faces',mesh_ctx.face,'vertices',mesh_ctx.vert),'FaceColor','b','EdgeColor','none','Parent',ha);

view(90,0); axis(ha,'image','off')
material shiny
lighting gouraud
hl = camlight('headlight');

% Headshape locations
%--------------------------------------------------------------------------
if ~isempty(headshape_polhemus)
    h_hsp   = plot3(headshape_polhemus(:,1),headshape_polhemus(:,2),headshape_polhemus(:,3),'dm');
    set(h_hsp,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','r');
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


function movelight(varargin)
    delete(hl)
    hl = camlight('headlight');
end

end