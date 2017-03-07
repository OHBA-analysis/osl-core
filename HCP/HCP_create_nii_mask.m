function HCP_create_nii_mask(sourcemodel,mask_out)
% Creates an nii mask of an HCP sourcemodel at the same resolution.
%
% HCP_create_nii_mask(sourcemodel,mask_out)
%
% sourcemodel   - filename of an HCP 3D sourcemodel
% mask_out      - filename of a nifti file to save to
%
% Adam Baker 2015

load(sourcemodel) % loads sourcemodel3d

[grid_x,grid_y,grid_z] = meshgrid(10*sourcemodel3d.cfg.grid.template.xgrid, ...
                                  10*sourcemodel3d.cfg.grid.template.ygrid, ...
                                  10*sourcemodel3d.cfg.grid.template.zgrid);
                              

grid_x = round(grid_x);                              
grid_y = round(grid_y);                              
grid_z = round(grid_z);                              
                              
xstep = mean(diff(sourcemodel3d.cfg.grid.template.xgrid*10));
ystep = mean(diff(sourcemodel3d.cfg.grid.template.ygrid*10));
zstep = mean(diff(sourcemodel3d.cfg.grid.template.zgrid*10));
     

coords_inside = round(10*sourcemodel3d.cfg.grid.template.pos(sourcemodel3d.cfg.grid.template.inside,:));

vol = zeros(size(grid_x));
for i = 1:numel(vol)
    v = sum(all(bsxfun(@minus,coords_inside,[grid_x(i),grid_y(i),grid_z(i)])==0,2));
    vol(i) = v;
end

vol = permute(vol,[2 1 3 4]);     
     
save_avw(vol,mask_out,'f',[xstep ystep zstep 1]);

xform = [xstep  0       0       round(10*min(sourcemodel3d.cfg.grid.template.xgrid)),...
         0      ystep   0       round(10*min(sourcemodel3d.cfg.grid.template.ygrid)),...
         0      0       zstep   round(10*min(sourcemodel3d.cfg.grid.template.zgrid)),...
         0      0       0       1];

runcmd(['fslorient -setsformcode 0 ' mask_out])
runcmd(['fslorient -setqformcode 2 ' mask_out])
runcmd(['fslorient -setsform ' num2str(xform) ' ' mask_out])
runcmd(['fslorient -setqform ' num2str(xform) ' ' mask_out])

end