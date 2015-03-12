function [dat4]=matrix2vols(dat2,mask)
%takes 2d matrix (space x time) and 3d mask (xyz)
%and returns 4d data (xyz x time)

%mask=reshape(mask,prod(size(mask)),1)'>0;
%dat4=zeros(prod(size(mask)),size(dat2,2));
%dat4(:,mask)=dat2;
%dat4=reshape(dat4,size(mask,1),size(mask,2),size(mask,3),size(dat2,2));

mask2=find(reshape(mask,prod(size(mask)),1)>0);
dat4=zeros(prod(size(mask)),size(dat2,2));
try
    dat4(mask2,:)=dat2;
catch
    error('mask and nvoxels in the data incompatible');
end;

dat4=reshape(dat4,size(mask,1),size(mask,2),size(mask,3),size(dat2,2));