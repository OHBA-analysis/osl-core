function [ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname)

% [ mni_coords xform ] = osl_mnimask2mnicoords(mask)
%
% converts an MNI standard brain mask, mask, into a list of mni_coords.

[ mni_res ] = get_nii_spatial_res( mask_fname );

mask=read_avw(mask_fname);

xstart=[90 -126 -72];
xform=eye(4)*1;xform(1,1)=-mni_res(1);xform(2,2)=mni_res(2);xform(3,3)=mni_res(3);
xform(1:3,4)=xstart;
mni_coord_vol=zeros([size(mask),3]);

for x=1:size(mask,1),
for y=1:size(mask,2),
for z=1:size(mask,3),

    if(mask(x,y,z)>0)
        tmp=round(xform*[([x-1 y-1 z-1]) 1]');
        mni_coord_vol(x,y,z,:)=[tmp(1:3)];
    end;

end
end
end

mni_coords=vols2matrix(mni_coord_vol,mask);
