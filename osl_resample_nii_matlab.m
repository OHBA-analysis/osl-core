function osl_resample_nii_matlab(input_nii,output_fname,output_mask,varargin)
    % Resample a nii file onto a different spatial grid using Matlab interpolation
    %
    % INPUTS
    %   input_nii : file name with nii volume to interpolate. Can also be standard MNI spatial resolution (e.g. 8
    %                If the file doesn't exist, try appending the osl standard mask folder
    %   output_mask: file name of template that the input should be interpolated onto. Can also be standard MNI spatial resolution (e.g. 8
    %                If the file doesn't exist, try appending the osl standard mask folder
    %   output_fname : save the result to this file
    %   key-value pairs
    %       - interptype : passed to flirt e.g. trilinear
    %       - enforce_mask : threshold the interpolated input with the output mask
    %       - force_positive : Assign 0 to any negative values (in case they are introduced by the resampling process)
    % 
    % EXAMPLE
    %   osl_resample_nii('in.nii.gz',8,'out.nii.gz','interptype','nearestneighbour','dilate',true)
    %
    % Romesh Abeysuriya 2017
    
    arg = inputParser;
    arg.addParameter('interptype','nearest'); % Anything supported by griddedInterpolant e.g. linear, cubic
    arg.addParameter('enforce_mask',true); % assign 0 to any values where the output mask is zero
    arg.addParameter('force_positive',false); % assign 0 to any negative values. Helpful if resampling a positive matrix where extrapolation could potentially produce negative values
    arg.parse(varargin{:});

    if isnumeric(input_nii)
        input_nii = fullfile(osldir,'std_masks',sprintf('MNI152_T1_%dmm_brain.nii.gz',input_nii));
    end

    if ~exist(input_nii)
        input_nii = fullfile(osldir,'std_masks',input_nii);
    end

    if isnumeric(output_mask)
        output_mask = fullfile(osldir,'std_masks',sprintf('MNI152_T1_%dmm_brain.nii.gz',output_mask));
    end

    if ~exist(output_mask)
        output_mask = fullfile(osldir,'std_masks',output_mask);
    end

    assert(exist(input_nii)~=0,'Input file not found')
    assert(exist(output_mask)~=0,'Output mask file not found')

    [input_x,input_y,input_z,input_vol,input_xform,input_step] = load_nii_vol(input_nii,true);
    [output_x,output_y,output_z,output_mask_vol,output_xform,output_step] = load_nii_vol(output_mask,false);

    output_vol = zeros([size(output_x),size(input_vol,4)]);

    for j = 1:size(input_vol,4)
        F = griddedInterpolant(input_x,input_y,input_z,input_vol(:,:,:,j),arg.Results.interptype);
        output_vol(:,:,:,j) = F(output_x,output_y,output_z);
    end

    if arg.Results.enforce_mask
        output_vol = bsxfun(@times,output_vol,output_mask_vol~=0);
    end

    if arg.Results.force_positive
        output_vol(output_vol < 0) = 0;
    end

    osl_save_nii(output_vol,output_step,output_xform,output_fname);

function [xg,yg,zg,vol,xform,step] = load_nii_vol(fname,flip)
    % Read the input nii file, and construct nd arrays for the spatial coordinates

    if nargin < 2 || isempty(flip) 
        flip = false;
    end
        
    [vol,input_res,xform] = osl_load_nii(fname);
    step = input_res.*diag(sign(xform(1:3,1:3)))';

    x = xform(1,4):step(1):(xform(1,4)+step(1)*(size(vol,1)-1));
    y = xform(2,4):step(2):(xform(2,4)+step(2)*(size(vol,2)-1));
    z = xform(3,4):step(3):(xform(3,4)+step(3)*(size(vol,3)-1));

    [xg,yg,zg] = ndgrid(x,y,z);

    if flip
        % Flip the orientation of dimensions so that the data is increasing 
        % in each dimension - in cases where xform is negative (like the OSL standard masks)
        % this flipping is required for griddedInterpolant to work
        if step(1) < 0
            xg = xg(end:-1:1,:,:);
            yg = yg(end:-1:1,:,:);
            zg = zg(end:-1:1,:,:);
            vol = vol(end:-1:1,:,:,:);
            xform(1,1) = -xform(1,1);
            xform(1,4) = x(end);
        end

        if step(2) < 0
            xg = xg(:,end:-1:1,:);
            yg = yg(:,end:-1:1,:);
            zg = zg(:,end:-1:1,:);
            vol = vol(:,end:-1:1,:,:);
            xform(2,2) = -xform(2,2);
            xform(1,4) = y(end);

        end

        if step(3) < 0
            xg = xg(:,:,end:-1:1,:);
            yg = yg(:,:,end:-1:1,:);
            zg = zg(:,:,end:-1:1,:);
            vol = vol(:,:,end:-1:1,:);
            xform(3,3) = -xform(3,3);
            xform(1,4) = z(end);
        end
        
    end
    

