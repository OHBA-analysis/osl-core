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

% Make the head upright
%Currently DOES NOT rotate the scatter objects
hc = get(ha,'Children');
local_rotate(ha,hc,[0 1 0],-35)

end

function local_rotate(ax,h,azel,alpha)
%% Patched version of 'rotate' from Matlab

% Determine the origin (center of plot box).
origin = sum([get(ax,'xlim')' get(ax,'ylim')' get(ax,'zlim')'])/2;

% find unit vector for axis of rotation
if numel(azel) == 2 % theta, phi
    theta = pi*azel(1)/180;
    phi = pi*azel(2)/180;
    u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
elseif numel(azel) == 3 % direction vector
    u = azel(:)/norm(azel);
end

alph = alpha*pi/180;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

for i=1:numel(h),
  t = get(h(i),'type');
  skip = 0;
  
  if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'patch') || strcmp(t,'scatter')
    
    % If patch, rotate vertices  
    if strcmp(t,'patch')
       verts = get(h(i),'Vertices');
       x = verts(:,1); y = verts(:,2); 
       if size(verts,2)>2
          z = verts(:,3);
       else
          z = [];
       end
       
    % If surface or line, rotate {x,y,z}data   
    else
       x = get(h(i),'xdata');
       y = get(h(i),'ydata');
       z = get(h(i),'zdata');
    end
    
    if isempty(z)
       z = -origin(3)*ones(size(y));
    end
    [m,n] = size(z);
    if numel(x) < m*n
      [x,y] = meshgrid(x,y);
    end
  elseif strcmp(t,'text')
    p = get(h(i),'position');
    x = p(1); y = p(2); z = p(3);
  elseif strcmp(t,'image')
    x = get(h(i),'xdata');
    y = get(h(i),'ydata');
    z = zeros(size(x));
  else
    skip = 1;
  end
  
  if ~skip,
    [m,n] = size(x);
    newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
    newxyz = newxyz*rot;
    newx = origin(1) + reshape(newxyz(:,1),m,n);
    newy = origin(2) + reshape(newxyz(:,2),m,n);
    newz = origin(3) + reshape(newxyz(:,3),m,n);

    if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'scatter')
      set(h(i),'xdata',newx,'ydata',newy,'zdata',newz);
    elseif strcmp(t,'patch')
      set(h(i),'Vertices',[newx,newy,newz]); 
    elseif strcmp(t,'text')
      set(h(i),'position',[newx newy newz])
    elseif strcmp(t,'image')
      set(h(i),'xdata',newx,'ydata',newy)
    end
  end
end
end


