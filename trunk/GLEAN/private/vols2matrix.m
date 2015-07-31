function [dat2]=vols2matrix(dat4,mask);
%[dat2,ind]=vols2matrix(dat4,mask);
%takes a 4d volume and mask and returns the 2d matrix 
%(space x time) 

mask=reshape(mask,prod(size(mask)),1)'>0;
dat2=reshape(dat4,prod(size(mask)),size(dat4,4))';
dat2=dat2(:,mask)';
