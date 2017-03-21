function HCP_mni2hcp(sourcemodel,fname_in,fname_out)
% Resample an MNI mask to the HCP grid, e.g. to create standard masks or
% parcellations based on MNI space templates.
%
% HCP_mni2hcp(sourcemodel,fname_in,fname_out)
%
% sourcemodel   - filename of an HCP 3D sourcemodel
% fname_in      - MNI space template or mask
% fname_out     - HCP grid template or mask
%
% Adam Baker 2015

load(sourcemodel) % loads sourcemodel3d

[vol_orig,dims,scales] = read_avw(fname_in);

%load(sourcemodel); % loads sourcemodel3d

mni_grid_x =   90 - (0:dims(1)-1)*scales(1);
mni_grid_y = -126 + (0:dims(2)-1)*scales(2);
mni_grid_z =  -72 + (0:dims(3)-1)*scales(3);

[grid_x,grid_y,grid_z] = meshgrid(10*sourcemodel3d.cfg.grid.template.xgrid, ...
                                  10*sourcemodel3d.cfg.grid.template.ygrid, ...
                                  10*sourcemodel3d.cfg.grid.template.zgrid);

vol_new = zeros([size(grid_x),size(vol_orig,4)]);

for p = 1:size(vol_orig,4)                              
    vol_new(:,:,:,p) = interp3(mni_grid_x, ...
                            mni_grid_y, ...
                            mni_grid_z, ...
                            permute(vol_orig(:,:,:,p),[2 1 3]), ...
                            grid_x, ...
                            grid_y, ...
                            grid_z, ...
                            'nearest');
end

vol_new = permute(vol_new,[2 1 3 4]);
     
xstep = mean(diff(sourcemodel3d.cfg.grid.template.xgrid*10));
ystep = mean(diff(sourcemodel3d.cfg.grid.template.ygrid*10));
zstep = mean(diff(sourcemodel3d.cfg.grid.template.zgrid*10));

save_avw(vol_new,fname_out,'f',[xstep ystep zstep 1]);

xform = [xstep  0       0       round(10*min(sourcemodel3d.cfg.grid.template.xgrid)),...
         0      ystep   0       round(10*min(sourcemodel3d.cfg.grid.template.ygrid)),...
         0      0       zstep   round(10*min(sourcemodel3d.cfg.grid.template.zgrid)),...
         0      0       0       1];
     
    
runcmd(['fslorient -setsformcode 0 ' fname_out])
runcmd(['fslorient -setqformcode 2 ' fname_out])
runcmd(['fslorient -setsform ' num2str(xform) ' ' fname_out])
runcmd(['fslorient -setqform ' num2str(xform) ' ' fname_out])


%runcmd(['fslview ' fname_out ' &'])
%runcmd(['fslview ' fname_in ' &'])
end
