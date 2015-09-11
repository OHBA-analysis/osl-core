function vec_out = HCP_mni2hcp(sourcemodel3d,fname_in,fname_out)
%HCP_mni2hcp resample mni nifti files onto the hcp grid
%
% vol_new = HCP_MNI2HCP(SOURCEMODEL, FNAMEIN) reads in nifti FNAMEIN and
%   uses hcp SOURCEMODEL to resample onto the hcp grid, selecting the 
%   internal parts.  
% 
% [] = HCP_MNI2HCP(SOURCEMODEL, FNAMEIN, FNAMEOUT) saves the output as a
%   nifti file with name FNAMEOUT

[vol_orig,dims,scales] = read_avw(fname_in);

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

% resampled volume in 3d hcp space
vol_new = permute(vol_new,[2 1 3 4]);
     
% turn into vector format
% fill in volume with source variance at correct locations
vec_out = zeros(length(sourcemodel3d.cfg.grid.template.inside),size(vol_orig,4));
for vox = sourcemodel3d.cfg.grid.template.inside;
   fillInd = (vox==sourcemodel3d.cfg.grid.template.inside);
   vec_out(fillInd,:) = ...
        vol_new(sourcemodel3d.cfg.grid.template.pos(vox,1) == sourcemodel3d.cfg.grid.template.xgrid,...
                sourcemodel3d.cfg.grid.template.pos(vox,2) == sourcemodel3d.cfg.grid.template.ygrid,...
                sourcemodel3d.cfg.grid.template.pos(vox,3) == sourcemodel3d.cfg.grid.template.zgrid,:);
end%for

%%% FUDGE FACTOR
vec_out(isnan(vec_out)) = 0;

% output as nifti
if nargin >= 3,
    xstep = mean(diff(sourcemodel3d.cfg.grid.template.xgrid*10));
    ystep = mean(diff(sourcemodel3d.cfg.grid.template.ygrid*10));
    zstep = mean(diff(sourcemodel3d.cfg.grid.template.zgrid*10));
    
    save_avw(vol_new,fname_out,'f',[xstep ystep zstep 1]);
    
    xform = [xstep  0       0       round(10*min(sourcemodel3d.cfg.grid.template.xgrid)),...
             0      ystep   0       round(10*min(sourcemodel3d.cfg.grid.template.ygrid)),...
             0      0       zstep   round(10*min(sourcemodel3d.cfg.grid.template.zgrid)),...
             0      0       0       1];
    
    
    dos(['fslorient -setsformcode 0 ' fname_out])
    dos(['fslorient -setqformcode 2 ' fname_out])
    dos(['fslorient -setsform ' num2str(xform) ' ' fname_out])
    dos(['fslorient -setqform ' num2str(xform) ' ' fname_out])
end
end

