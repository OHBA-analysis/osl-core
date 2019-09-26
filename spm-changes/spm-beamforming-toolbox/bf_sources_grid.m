function res = bf_sources_grid(BF, S)
% Generate beamforming grid
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_sources_grid.m 117 2014-07-17 13:11:48Z litvak.vladimir@gmail.com $

%--------------------------------------------------------------------------
if nargin == 0 
    resolution = cfg_entry;
    resolution.tag = 'resolution';
    resolution.name = 'Grid resolution';
    resolution.strtype = 'n';
    resolution.num = [1 1];
    resolution.val = {5};
    resolution.help = {'Select the resolution of the grid (in mm)'};
    
    space = cfg_menu;
    space.tag = 'space';
    space.name = 'Coordinate system';
    space.help = {'Select the coordinate system in which the grid should be generated'};
    space.labels = {'MNI template', 'MNI-aligned', 'Head', 'Native'};
    space.values = {'MNI template', 'MNI-aligned', 'Head', 'Native'};
    space.val = {'MNI template'};

    grid = cfg_branch;
    grid.tag = 'grid';
    grid.name = 'Grid';
    grid.val = {resolution, space};
    
    res = grid;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

iskull = export(gifti(BF.data.mesh.tess_iskull), 'ft');

M1 = BF.data.transforms.toNative;

switch S.space
    case 'MNI template'
        M1 = BF.data.transforms.toMNI/M1;
        M2 = inv(BF.data.transforms.toMNI);
    case 'MNI-aligned'
        M1 = BF.data.transforms.toMNI_aligned/M1;
        M2 = inv(BF.data.transforms.toMNI_aligned);
    case 'Head'
        M1 = BF.data.transforms.toHead/M1;
        M2 = inv(BF.data.transforms.toHead);
    case 'Native'
        M2 = inv(M1);
        M1 = eye(4);
end

iskull = ft_convert_units(ft_transform_geometry(M1, iskull));

mn = min(iskull.pnt);
mx = max(iskull.pnt);

resolution = S.resolution;

if isequal(iskull.unit, 'm')
    resolution = 1e-3*resolution;
end

% If zero is inside the brain, make sure grid points fall on multiples of
% resolution to ease simulating data from points on the grid
if mn(1)<0 && mx(1)>0
    grid.xgrid = [fliplr(0:-resolution:mn(1)) resolution:resolution:mx(1)];
else
    grid.xgrid = mn(1):resolution:mx(1);
end

if mn(2)<0 && mx(2)>0
    grid.ygrid = [fliplr(0:-resolution:mn(2)) resolution:resolution:mx(2)];
else
    grid.ygrid = mn(2):resolution:mx(2);
end

if mn(3)<0 && mx(3)>0
    grid.zgrid = [fliplr(0:-resolution:mn(3)) resolution:resolution:mx(3)];
else
    grid.zgrid = mn(3):resolution:mx(3);
end

grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
[X, Y, Z]  = ndgrid(grid.xgrid, grid.ygrid, grid.zgrid);

pos   = [X(:) Y(:) Z(:)];

inside = ft_inside_vol(pos, struct('bnd', iskull));

pos    = spm_eeg_inv_transform_points(M2, pos);

grid.allpos  = pos;
grid.inside  = find(inside);
grid.outside = find(~inside);

grid.pos     = pos(inside, :);

res = grid;