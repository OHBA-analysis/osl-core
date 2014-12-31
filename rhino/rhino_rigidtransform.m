function [x,Mr] = rhino_rigidtransform(x,rx,ry,rz,tx,ty,tz)

Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
T  = [tx ty tz];
Mr = [Rx*Ry*Rz T'; 0 0 0 1];

x  = Mr*[x; ones(1,size(x,2))];
x  = x(1:3,:);

end