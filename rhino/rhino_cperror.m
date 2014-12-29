function dists = rhino_cperror(q,p,method)
% Compute distances between points p and matching points in point cloud q
% E.g. q = MRI mesh, p = Polhemus points
%
% dist = RHINO_CPERROR(q,p)
%
% AB 2014

if size(p,2)==3
    p = p';
end

if size(q,2)==3
    q = q';
end

if nargin < 3
    method = 'surface';
end

switch method
    case 'point'
        dists = point_to_point_distance(q,p);
        
    case 'surface'
        dists = point_to_surface_distance(q,p);
end

dists = dists(:);

end


function dists = point_to_surface_distance(q,p)
% Get closest points in q to p to define each triangular surface patch:
neighbours = knnsearch(q',p','k',3);

% Loop through each matching point and compute closest distance from point
% in p to each triangular patch
dists = zeros(size(neighbours,1),1);
for i = 1:length(dists)
    dists(i) = pointTriangleDistance(q(:,neighbours(i,:))',p(:,i)');
end

end




function dists = point_to_point_distance(q,p)
m = size(p,2);
n = size(q,2);
dists = zeros(1,m);
for ki=1:m
    d=zeros(1,n);
    for ti=1:3
        d=d+(q(ti,:)-p(ti,ki)).^2;
    end
    dists(ki) = min(d);
end

dists = sqrt(dists);

end

