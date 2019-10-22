function [ fname_out, lost_coords_indices ] = osl_mnicoords2mnimask( mni_coords_in, gridstep, fname, resamp_gridstep, roi_radius )

% [ fname_out ] = osl_mnicoords2mnimask.m( mni_coords_in, gridstep, fname, resamp_gridstep )
%
% creates a mask using the passed in mni coordinates. Create a mask with
% voxels on a grid with spacings of gridstep, and then downsamples this to
% resamp_gridstep
%
% if roi_radius=-1 then finds single nearest point in MNI152_T1_' num2str(gridstep) 'mm_brain mask 
% lost_coords_indices only calculated if roi_radius=-1

OSLDIR = getenv('OSLDIR');

if(nargin<4)
    resamp_gridstep=2;
end

if(nargin<5)
    roi_radius=-1;
end

mask_fname=[OSLDIR '/std_masks/MNI152_T1_' num2str(gridstep) 'mm_brain.nii.gz'];
[mask,~,xformmni] = nii.load(mask_fname);
newmask=zeros(size(mask));
newmask=vols2matrix(newmask,mask);

lost_coords_indices=zeros(size(mni_coords_in,1),1);
dist_store=zeros(size(newmask));

[ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname);

for mm=1:size(mni_coords_in,1),    
    dists=(sqrt(sum((mni_coords-repmat(mni_coords_in(mm,:),length(mni_coords),1)).^2,2)));
    
    if(roi_radius<0)
        [dist,seed_index]=min(dists); 
    
        %disp(['Nearest grid point at MNI coordinate '
        %num2str(mni_coords(seed_index,:))]);

        if newmask(seed_index)==0 && dist < gridstep
            
            newmask(seed_index)=1;
            dist_store(seed_index)=dist;
            
        else
           
            lost_coords_indices(mm)=1;
            
        end
    else
        newmask(find(dists<roi_radius))=1;
    end
    
end

newmask2=matrix2vols(newmask,mask);

nii.save(newmask2,gridstep,xformmni,[fname '_' num2str(gridstep) 'mm']); 

if(resamp_gridstep~=gridstep)
    fname_out = nii.resample([fname '_' num2str(gridstep) 'mm'],[fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'interptype','cubic','enforce_mask',true);
else
    fname_out=[fname '_' num2str(gridstep) 'mm'];
end