function fsleyes(fnames,thresholds,colour_maps,anatomical,filename_for_fsleyes)
    % Wrapper to call fsleyes
    %
    % INPUTS
    % - fnames - A file name or cell array of file names (mandatory)
    % - thresholds - Array with two elements, or cell array of such elements
    %   the same size as fnames with colour limits (default: [] will do no thresholding)
    % - colour_maps - String, or cell array of strings, with name of valid
    %   colour map (default: 'red-yellow', specify 'greyscale' for grey scale! )
    % - anatomical - File name of base image to use (default: automatic
    % std, specify 'none' if you want no anatomical file loaded), 
    %
    % e.g.: To show 2 images with no thresholding, with greyscale colormaps and no base anatomical image:
    %       fsleyes({niftii1, niftii2},[],'greyscale','none')
    %
    % Romesh Abeysuriya 2017

    
    if nargin < 5 || isempty(filename_for_fsleyes)
        filename_for_fsleyes='fsleyes';
    
        filename_for_fsleyes='/Applications/FSLeyes.app/Contents/MacOS/fsleyes';    
    end
    
    assert(nargin > 0,'You must specify at least one image file to display');
    assert(ischar(fnames) || iscell(fnames),'Input must be either a file name or a cell array of file names');

    OSLDIR = getenv('OSLDIR');

    available_maps = {'blue-lightblue','blue','cool','copper','cortical','green','greyscale','hot','pink','random','red-yellow','red','render1','render1t','render2','render2t','render3','subcortical','yellow'};

    if nargin < 4 || isempty(anatomical)
        anatomical='';
    end

    if nargin < 3 || isempty(colour_maps)
        colour_maps = 'red-yellow';
    end

    if nargin < 2 ||isempty(thresholds)
         thresholds = [];
    end

    % Deal with expanding all inputs into cell arrays and extending thresholds and colour maps if required
    if ischar(fnames)
        fnames = {fnames};
    end

    % Check all files exist
    for j = 1:length(fnames)
        if ~osl_util.isfile(fnames{j})
            error(sprintf('File not found: %s',fnames{j}));
        end
    end

    if isnumeric(thresholds)
        thresholds = {thresholds};
    end

    if ischar(colour_maps)
        colour_maps = {colour_maps};
    end

    assert(length(thresholds)==1 || length(thresholds) == length(fnames),'Number of threshold pairs must be scalar or the same as the number of image files');
    assert(length(colour_maps)==1 || length(colour_maps) == length(fnames),'Number of threshold pairs must be scalar or the same as the number of image files');

    if length(thresholds) ~= length(fnames) 
        thresholds(1:length(fnames)) = thresholds(1);
    end

    if length(colour_maps) ~= length(fnames) 
        colour_maps(1:length(fnames)) = colour_maps(1);
    end

    % Construct the command
    file_string = '';
    res = nan(length(fnames),1);
    for j = 1:length(fnames)
        tmp = nii.get_spatial_res(fnames{j});
        res(j) = tmp(1);

        assert(any(ismember(colour_maps{j},available_maps)),'Colour map %s not recognized',colour_maps{j});
        file_string = sprintf('%s %s -cm %s',file_string,fnames{j},colour_maps{j});
        if ~isempty(thresholds{j})
            assert(length(thresholds{j})==2,'Display range must be empty or two numbers')
            file_string = sprintf('%s -dr %.2f %.2f',file_string,thresholds{j}(1),thresholds{j}(2));
        end
    end

    assert(numel(unique(res)) == 1,'Spatial maps must be of equal sizes')
    
    if isempty(anatomical) && ~strcmp(anatomical,'none')
        anatomical = fullfile(osldir,'std_masks',sprintf('MNI152_T1_%dmm_brain',res(1))); 
    else
        anatomical='';
    end

    % Current version of fsleyes returns 0 even if an error occurred. So this command
    % below will fail silently. Hopefully this will be fixed upstream later on
    runcmd('%s %s %s &',filename_for_fsleyes, anatomical,file_string);
     
