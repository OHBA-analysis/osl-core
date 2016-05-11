function [img,dims,scales,bpp,endian] = ra(fname)
% [img, dims,scales,bpp,endian] = RA(fname)
%
% calls read_avw
%
%  See also: READ_AVW SAVE_AVW

[img,dims,scales,bpp,endian] = read_avw(fname);