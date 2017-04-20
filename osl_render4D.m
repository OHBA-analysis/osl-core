function osl_render4D(nii,varargin)
    % Creates a surface rendering of a 4D nifti file and saves as dense time series
    % (.dtseries.nii) CIFTI file using HCP workbench
    %
    % INPUTS
    %   osl_render4d(nii,varargin)
    %
    % where nii is the name of the nii file with a 4D matrix to render (include the extension)
    %
    % where varargin corresponds to key-value pairs for
    % 'savedir' : location to save output files (default=current directory)
    % 'interptype' : 'trilinear' or for nearest-neighbour, 'enclosing' (default='trilinear')
    % 'visualise' : show using workbench (default=true)
    %
    % EXAMPLE USAGE
    %   osl_render_4d('power_map.nii.gz','interptype','enclosing','visualise',false)
    % 
    % Romesh Abeysuriya 2017
    % Adam Baker 2013

    % Input nii may or may not have extension
    if ~exist(nii) || isempty(strfind(nii,'.nii'))
      error('input should be nii or nii.gz file');
    end
    
    [inpath,infile] = fileparts(nii);

    arg = inputParser;
    arg.addParameter('savedir',inpath); 
    arg.addParameter('interptype','trilinear'); 
    arg.addParameter('visualise',true); 
    arg.parse(varargin{:});

    if ~isempty(arg.Results.savedir) && ~isdir(arg.Results.savedir)
        mkdir(arg.Results.savedir);
    end

    infile = strrep(infile,'.gz','');
    infile = strrep(infile,'.nii','');
    outfile = fullfile(arg.Results.savedir,infile);


    % Load surfaces to map to
    surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
    surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');
    surf_right_inf = fullfile(osldir,'std_masks','ParcellationPilot.R.inflated.32k_fs_LR.surf.gii');
    surf_left_inf = fullfile(osldir,'std_masks','ParcellationPilot.L.inflated.32k_fs_LR.surf.gii');
    surf_right_vinf = fullfile(osldir,'std_masks','ParcellationPilot.R.very_inflated.32k_fs_LR.surf.gii');
    surf_left_vinf = fullfile(osldir,'std_masks','ParcellationPilot.L.very_inflated.32k_fs_LR.surf.gii');

    output_right    = [outfile '_right.func.gii'];
    output_left     = [outfile '_left.func.gii'];

    % Map volume to surface
    runcmd('wb_command -volume-to-surface-mapping %s %s %s -%s',nii,surf_right,output_right,arg.Results.interptype)
    runcmd('wb_command -volume-to-surface-mapping %s %s %s -%s',nii,surf_left,output_left,arg.Results.interptype)

    % Save as dtseries 
    cifti_right = strrep(output_right,'.func.gii','.dtseries.nii');
    cifti_left  = strrep(output_left, '.func.gii','.dtseries.nii');

    runcmd('wb_command -cifti-create-dense-timeseries %s -right-metric %s',cifti_right,output_right)
    runcmd('wb_command -cifti-create-dense-timeseries %s -left-metric %s',cifti_left,output_left)

    % View in workbench
    if arg.Results.visualise
      runcmd('wb_view %s %s %s %s %s %s %s %s &',surf_left,surf_right,surf_left_inf,surf_right_inf,surf_left_vinf,surf_right_vinf,cifti_left,cifti_right);
    end
