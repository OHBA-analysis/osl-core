function [bin,vecx,vecy] = hist2d(mat,x,y);

% bin = hist2d(mat);
% 2-D histogram with 10*10 bins
% mat should be N*2, where first col is x and second is y
%
% bin = hist2d(mat,n);
% 2-D histogram with n*n bins
%
% bin = hist2d(mat,nx,ny);
% 2-D histogram with nx*ny bins
%
% bin = hist2d(mat,x);
% 2-D histogram with bins specified by x in both directions 
% (x must be regularly spaced)
%
% bin = hist2d(mat,x,y);
% 2-D histogram with bins specified by x and y 
% (x and y must be regularly spaced)
%
% Copyright(C) 2001
%  T.Behrens, M.Jenkinson
%

if(size(mat,2) ~=2), mat=mat'; end;

if(nargin == 1),
  x=10;
end;

if(length(x) == 1)
  nx=x;
  
  if(nargin <= 2),
    ny=x;
  else,
    ny=y;
  end;

  maxx=max(mat(:,1));minx=min(mat(:,1));rx=maxx-minx;
  maxy=max(mat(:,2));miny=min(mat(:,2));ry=maxy-miny;
else,
  nx = length(x);	
  maxx=max(x);minx=min(x);rx=maxx-minx;
  if(nargin == 2),
    ny=length(x);
    maxy=max(x);miny=min(x);ry=maxy-miny;
  else,
    ny=length(y);	
    maxy=max(y);miny=min(y);ry=maxy-miny;
  end;
end;

bin = zeros(nx,ny);

for i=1:length(mat);
  xbin=ceil(nx*(mat(i,1)-minx)/rx);
  ybin=ceil(ny*(mat(i,2)-miny)/ry);
  if(~(xbin<1 | ybin<1 | xbin>nx | ybin>ny))
    bin(xbin,ybin)=bin(xbin,ybin)+1;
  end;
end
vecy=[miny:ry/(ny-1):maxy];
vecx=[minx:rx/(nx-1):maxx];
%keyboard
if(nargout==0)
  imagesc(vecy,vecx,bin(1:size(bin,1),1:size(bin,2)));
  colorbar;
  axis xy;axis image;
end










