function [ fname_out ] = osl_mnicoords2mnimask( mni_coords_in, gridstep, fname, resamp_gridstep, roi_radius )

% [ fname_out ] = osl_mnicoords2mnimask.m( mni_coords_in, gridstep, fname, resamp_gridstep )
%
% creates a mask using the passed in mni coordinates. Create a mask with
% voxels on a grid with spacings of gridstep, and then downsamples this to
% resamp_gridstep

OSLDIR = getenv('OSLDIR');

if(nargin<4)
    resamp_gridstep=2;
end;

if(nargin<5)
    roi_radius=-1;
end;

mask_fname=[OSLDIR '/std_masks/MNI152_T1_' num2str(gridstep) 'mm_brain'];
mask=read_avw(mask_fname);
newmask=zeros(size(mask));
newmask=vols2matrix(newmask,mask);

[ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname);

for mm=1:size(mni_coords_in,1),    
    dists=(sqrt(sum((mni_coords-repmat(mni_coords_in(mm,:),length(mni_coords),1)).^2,2)));
    
    if(roi_radius<0)
        [dist,seed_index]=min(dists); 
    
        %disp(['Nearest grid point at MNI coordinate '
        %num2str(mni_coords(seed_index,:))]);

        newmask(seed_index)=1;
    else
        newmask(find(dists<roi_radius))=1;
    end;
end;

newmask2=matrix2vols(newmask,mask);

save_avw(newmask2,[fname '_' num2str(gridstep) 'mm'], 'f', [gridstep,gridstep,gridstep,1]); 

if(resamp_gridstep~=gridstep)
    osl_resample_nii([fname '_' num2str(gridstep) 'mm'],[fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[OSLDIR '/std_masks/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
    fname_out=[fname '_' num2str(resamp_gridstep) 'mm'];
else
    fname_out=[fname '_' num2str(gridstep) 'mm'];
end;