function fname_out = nii_quicksave(mat,fname,options_or_input_spat_res,output_spat_res,interp)

% fname_out = nii_quicksave(mat,fname,options)
% OR
% fname_out = nii_quicksave(mat,fname,input_spat_res,output_spat_res,interp) - interface provided for some
% backwards compatibility
% 
% Saves a niftii file for passed in data
%
% mat is a 2D matrix of data (nvoxels x ntpts)
%
% fname is the output filename
%
% options.output_spat_res is is output spatial resolution in mm
% options.interp is interpolation method to use for flirt resampling, {default 'trilinear'}
% options.tres is temporal res in secs {default 1s}
% options.mask_fname is niftii filename of mask used, if not passed in (or set to []) then a whole brain mask is assume
% 
% MWW 2015

global OSLDIR;

if nargin<3
    options_or_input_spat_res=struct();
end;

%%%%%%%%%%%%%%%
% backwards compatibility
if ~isstruct(options_or_input_spat_res)
    % old interface was:
    % nii_quicksave(mat,fname,input_spat_res,output_spat_res,interp)
    input_spat_res=options_or_input_spat_res;
    
    options=struct();
    if nargin>3
        options.output_spat_res=output_spat_res;
    end;
    if nargin>4
        options.interp=interp;
    end;
else
    options=options_or_input_spat_res;
end;

clear options_or_input_spat_res;

%%%%%%%%%%%%%%%
% parse options
try 
    options = ft_checkopt(options,'interp','char',{'trilinear','nearestneighbour','sinc','spline'});
catch
    options.interp='trilinear';
end;
interp=options.interp;

try 
    options = ft_checkopt(options,'tres','double');
catch
    options.tres=1;
end;
tres=options.tres;

try 
    options = ft_checkopt(options,'mask_fname','char');
catch
    options.mask_fname=[];
end;
mask_fname=options.mask_fname;

%%%%%%%%%%%%%%%
% establish input_spat_res and setup mask
if isempty(mask_fname)    
    % assume mask is wholebrain std space mask

    % establish input_spat_res 
    if ~exist('input_spat_res','var')
        input_spat_res=getmasksize(size(mat,1));
    else
        if input_spat_res~=getmasksize(size(mat,1))
            error('mat nvoxels and input_spat_res incompatible');
        end;
    end;
    
    % load in std brain mask
    mask_fname=[OSLDIR '/std_masks/MNI152_T1_' num2str(input_spat_res) 'mm_brain.nii.gz'];
    
    % for a sanity check:
    mask_spat_res = get_nii_spatial_res( mask_fname );
    mask_spat_res = mask_spat_res(1);
        
else
    % establish input_spat_res from mask hdr
    mask_spat_res = get_nii_spatial_res( mask_fname );
    mask_spat_res = mask_spat_res(1);
    input_spat_res = mask_spat_res;
end;

if input_spat_res~=mask_spat_res,
    error('input_spat_res~=mask_spat_res');
end;
    
% load in mask
stdbrain=read_avw(mask_fname); 

%%%%%%%%%%%%%%%
% set output spat res
try 
    options = ft_checkopt(options,'output_spat_res','double');
catch
    options.output_spat_res = input_spat_res;
end;
output_spat_res=options.output_spat_res;

%%%%%%%%%%%%%%%
% now do the actual work

% save nii
save_avw(matrix2vols(mat,stdbrain),fname,'f',[input_spat_res input_spat_res input_spat_res tres]);

% resample using flirt (this will also help make sure hdr info is set)
fname = strrep(fname,'.gz','');
fname = strrep(fname,'.nii','');
fname_rs=[fname '_ds' num2str(output_spat_res) 'mm'];
stdbrain = [OSLDIR '/std_masks/MNI152_T1_' num2str(output_spat_res) 'mm_brain.nii.gz'];
osl_resample_nii(fname, fname_rs, output_spat_res,interp,stdbrain);

% tidy up files
dos(['rm ' fname '.nii.gz']);
dos(['mv ' fname_rs '.nii.gz ' fname '.nii.gz']);

fname_out=[fname '.nii.gz'];

end
