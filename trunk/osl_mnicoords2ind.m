function [ ind ] = osl_mnicoords2ind(mni_coords, mni_res)

% [ ind ] = osl_mnicoords2ind(mni_coords, mni_res)
%
% converts an MNI coord into indexes into volume (starting from 0)

xstart=[90 -126 -72];
xform=eye(4)*1;xform(1,1)=-mni_res;xform(2,2)=mni_res;xform(3,3)=mni_res;
xform(1:3,4)=xstart;

ind=round(inv(xform)*[([mni_coords]) 1]');
ind=ind(1:3);
