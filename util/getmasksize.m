function s = getmasksize(n)
% Really dumb way of determining size of standard space mask knowing only
% the number of voxels in the grid (& vice versa)

if n == 676;
  s = 14;
elseif n == 14
  s = 676;
elseif n == 1065;
  s = 12;
elseif n == 12
  s = 1065;
elseif n == 1821;
  s = 10;
elseif n == 10
  s = 1821;
elseif n == 3559;
  s = 8;
elseif n == 8
  s = 3559;
elseif n == 8471
  s = 6;  
elseif n == 6
  s = 8471;    
elseif n == 228453;
  s = 2;
elseif n == 2
  s = 228453;
else
  error('unknown mask size')
end


end

