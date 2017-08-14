%% osl_save_nii_ica_maps.m
%  
%  This function produces nifti file maps of the spatial topographies of
%  ICA outputs.
%
%  Syntax: map_name = osl_save_nii_ica_maps(oil,type,order)
%
%  INPUTS: Required: oil   -  the OIL structure output by osl_run_oil.m
%          Optional: type  -  the type of map produced.
%                             'correlation' - default map for MEG ICA. A map
%                             of correlation values are produced between each
%                             IC and the concatenated input data.
%                             'scaled_covariance' - A map of of the
%                             covariance between each IC and the concatenated
%                             input data at each voxel is produced. Each map
%                             is scaled to have unit variance.
%                             'spatial' - Spatial maps from the spatial ICA
%                             decompostion.
%                             'variance' - gives average varaince map of beamformed data
%                     order - User-specified re-ordering of the ICs
%                             (usually based on the GLM statistics). order
%                             is a vector of length equal to
%                             oil.ica.num_ics.
%
%  OUTPUTS: map_name: the full path name of the nifit file. The files are
%                     always save in the OIL directory/ICA subdirectory
%                     [oil.name '/ica_dir/'].
%
%  HL 160113
% Version 1.0

%%
function map_name = oil_save_nii_ica_maps(oil,type,order)

OSLDIR = getenv('OSLDIR');

% parse inputs
if nargin<2,  type   =  'correlation'; try order=1:oil.ica.num_ics; catch, order=1:oil.num_ics;  end; end
if nargin<3,  try order=1:oil.ica.num_ics; catch, order=1:oil.num_ics;  end; end

% find save directory
try
save_dir = [oil.source_recon.dirname '/' oil.enveloping.name '/' oil.concat_subs.name '/' oil.ica.name];
if ~isdir(save_dir), mkdir(save_dir);end;
catch
[pat,nam]=fileparts(oil.ica_concat_path{1});
save_dir=[pat '/' oil.name];
if ~isdir(save_dir), mkdir(save_dir);end;
end

% get GM mask
if isfield(oil, 'use_gm_mask') && oil.use_gm_mask,
    try
        gm_mask = nii.quickread([OSLDIR '/std_masks/grey_matter/MNI_greymatter_priors_' num2str(oil.enveloping.gridstep) 'mm.nii.gz'],oil.enveloping.gridstep);
    catch
        warning('Unable to find Grey Matter Masks');
    end
end%if

switch type
    case {'spatial'}
        map_stem = 'ica_spatial_maps_';
        fname    = fullfile(save_dir, ...
                            [map_stem oil.concat_subs.results.concat_file]);
        map_name = nii.quicksave(oil.ica.results.sICs', fname, ...
                                 oil.enveloping.gridstep,2);
        
    case {'variance'}
        varMapName =  ['session', ...
                       num2str(oil.concat_subs.sessions_to_do(1)), ...
                       '_recon_noise_corrected_variance_map'];

        Nvoxels = numel(nii.quickread(fullfile(oil.source_recon.dirname, ...
                                               oil.enveloping.name, ...
                                               'variance_maps', ...
                                               varMapName), ...
                                      oil.source_recon.gridstep));
        varmaps = zeros(Nvoxels,1);

        for i = 1:numel(oil.concat_subs.sessions_to_do)
            singsub_varmaps_name = ['session', ...
                                    num2str(oil.concat_subs.sessions_to_do(i)), ...
                                    '_recon_noise_corrected_variance_map'];
            var_maps_singsub     = nii.quickread(singsub_varmaps_name, ...
                                                 oil.source_recon.gridstep);

            varmaps = varmaps + var_maps_singsub/std(var_maps_singsub);
        end

        varmaps = varmaps / numel(oil.concat_subs.sessions_to_do);
        fname   = fullfile(oil.source_recon.dirname, ...
                                               oil.enveloping.name, ...
                                               'variance_maps', ...
                                               'average_of_sessions_recon_noise_corrected_variance_map');
        map_name = nii.quicksave(varmaps, fname, ...
                                 oil.source_recon.gridstep, 2);
        
    case {'correlation'}
        
        map_stem='ica_correlation_maps_';
        if  isfield(oil,'concat_subs') && isfield(oil,'ica') && isfield(oil,'enveloping') && isfield(oil.ica.results,'tICs') && isfield(oil.concat_subs.results,'concat_file') && isfield(oil.enveloping,'gridstep')
            
            ica_concat_fname = fullfile(oil.source_recon.dirname, ...
                                        oil.enveloping.name, ...
                                        oil.concat_subs.name, ...
                                        oil.concat_subs.results.concat_file);

            try
                ica_concat = nii.quickread(ica_concat_fname, ...
                                           oil.enveloping.gridstep);
            catch % probably .mat file
                ica_concat = load(ica_concat_fname);
                ica_concat = ica_concat.ica_concat;
            end

            map = zeros(size(ica_concat,1),size(oil.ica.results.tICs,1));

            for i = 1:size(oil.ica.results.tICs,1)
                for v = 1:size(map,1)
                    map(v,i) = corr(oil.ica.results.tICs(order(i),:)', ica_concat(v,:)');
                end
            end

            fname    = [save_dir '/' map_stem oil.concat_subs.results.concat_file];
            map_name = nii.quicksave(map,fname,oil.enveloping.gridstep,2);

        elseif isfield(oil,'name') && isfield(oil,'ica_concat_path') && isfield(oil,'gridstep') && isfield(oil,'tICs')
            
            ica_concat_fname = oil.ica_concat_path{1};
            try
                ica_concat = nii.quickread(ica_concat_fname, ...
                                           oil.enveloping.gridstep);
            catch % probably .mat file
                ica_concat = load(ica_concat_fname);
                ica_concat = ica_concat.ica_concat;
            end

            map=zeros(size(ica_concat,1),size(oil.tICs,1));

            for i=1:size(oil.tICs,1)
                for v=1:size(map,1)
                    map(v,i)=corr(oil.tICs(order(i),:)',ica_concat(v,:)');
                end
            end

            [~, nam] = fileparts(oil.ica_concat_path{1});
            fname     = fullfile(save_dir , [map_stem nam]);

            if oil.use_gm_mask; map(gm_mask<0.3)=0; end

            map_name  = nii.quicksave(map,fname,oil.gridstep,2);
        else
            error('Unexpected input format')
        end
        
    case {'scaled_covariance'}
        
        map_stem='ica_covariance_maps_';
        if  isfield(oil,'concat_subs') && isfield(oil,'ica') && isfield(oil,'enveloping') && isfield(oil.ica.results,'tICs') && isfield(oil.concat_subs.results,'concat_file') && isfield(oil.enveloping,'gridstep')
            
            ica_concat_fname = fullfile(oil.source_recon.dirname, ...
                                        oil.enveloping.name, ...
                                        oil.concat_subs.name, ...
                                        oil.concat_subs.results.concat_file);
            try
                ica_concat = nii.quickread(ica_concat_fname, ...
                                           oil.enveloping.gridstep);
            catch % probably .mat file
                ica_concat = load(ica_concat_fname);
                ica_concat = ica_concat.ica_concat;
            end

            map = zeros(size(ica_concat,1),size(oil.ica.results.tICs,1));

            for i = 1:size(oil.ica.results.tICs,1)
                for v = 1:size(map,1)
                    map(v,i) = oil.ica.results.tICs(order(i),:) * demean(ica_concat(v,:)');
                end
                map(:,i)=map(:,i)/std(map(:,i));
            end

            fname    = [save_dir filesep map_stem oil.concat_subs.results.concat_file];
            map_name = nii.quicksave(map, fname, oil.enveloping.gridstep, 2);

        elseif isfield(oil,'name') && isfield(oil,'ica_concat_path') && isfield(oil,'gridstep') && isfield(oil,'tICs')
            
            ica_concat_fname = oil.ica_concat_path{1};
            try
                ica_concat = nii.quickread(ica_concat_fname, ...
                                           oil.enveloping.gridstep);
            catch % probably .mat file
                ica_concat = load(ica_concat_fname);
                ica_concat = ica_concat.ica_concat;
            end

            map = zeros(size(ica_concat,1), size(oil.tICs,1));

            for i = 1:size(oil.tICs,1)
                for v = 1:size(map,1)
                    map(v,i) = oil.tICs(order(i),:) * ica_concat(v,:)';
                end
                map(:,i) = map(:,i) / std(map(:,i));
            end
            
            [pat,nam] = fileparts(oil.ica_concat_path{1});
            fname     = fullfile(pat, [map_stem nam]);

            if oil.use_gm_mask; map(gm_mask<0.3)=0; end

            map_name  = nii.quicksave(map, fname, oil.gridstep, 2);
        else
            error('Unexpected input format')
        end
        
    case {'multiband_correlation'}
        
        map_stem = 'ica_correlation_maps_';
        if isfield(oil,'concat_subs') && isfield(oil,'ica') && isfield(oil,'enveloping') && isfield(oil.ica,'tICs') && isfield(oil.concat_subs,'concat_file') && isfield(oil.enveloping,'gridstep') && isfield(oil.ica,'tICs')
            for f = 1:size(oil.ica.ica_concat_path,1)
               
                ica_concat_fname = oil.ica.ica_concat_path{f};
                try
                    ica_concat = nii.quickread(ica_concat_fname, ...
                                               oil.enveloping.gridstep);
                catch % probably .mat file
                    ica_concat = load(ica_concat_fname);
                    ica_concat = ica_concat.ica_concat;
                end

                map = zeros(size(ica_concat,1), size(oil.ica.tICs,1));

                for i = 1:size(oil.ica.tICs,1)
                    for v = 1:size(map,1)
                        map(v,i) = corr(oil.ica.tICs(order(i),:)', ...
                                        ica_concat(v,:)');
                    end
                end

                [~,nam,~]   = fileparts(oil.ica.ica_concat_path{f});
                [~,nam,~]   = fileparts(nam);
                fname       = [save_dir '/' map_stem nam];

                if oil.ica.use_gm_mask; map(gm_mask<0.3)=0; end

                map_name{f} = nii.quicksave(map, fname, ...
                                            oil.enveloping.gridstep, 2);
            end
        else
            error('Unexpected input format')
        end
        
    case {'multiband_scaled_covariance'}
        
        map_stem = 'ica_covariance_maps_';
        if isfield(oil,'concat_subs') && isfield(oil,'ica') && isfield(oil,'enveloping') && isfield(oil.ica,'tICs') && isfield(oil.concat_subs,'concat_file') && isfield(oil.enveloping,'gridstep') && isfield(oil.ica,'tICs')
            for f = 1:size(oil.ica.ica_concat_path,1)
                
                ica_concat_fname = oil.ica.ica_concat_path{f};
                try
                    ica_concat = nii.quickread(ica_concat_fname, ...
                                               oil.enveloping.gridstep);
                catch % probably .mat file
                    ica_concat = load(ica_concat_fname);
                    ica_concat = ica_concat.ica_concat;
                end

                map = zeros(size(ica_concat,1), size(oil.ica.tICs,1));

                for i = 1:size(oil.ica.tICs,1)
                    for v = 1:size(map,1)
                        map(v,i) = oil.ica.tICs(order(i),:) * demean(ica_concat(v,:)');
                    end
                    map(:,i) = map(:,i) / std(map(:,i));
                end

                [~,nam,~]   = fileparts(oil.ica.ica_concat_path{f});
                [~,nam,~]   = fileparts(nam);
                fname       = [save_dir filesep map_stem nam];

                if oil.ica.use_gm_mask; map(gm_mask<0.3)=0; end

                map_name{f} = nii.quicksave(map, fname, oil.enveloping.gridstep, 2);
            end
        else
            error('Unexpected input format')
        end
        
        
        
    otherwise
        error('Unexpected ica type. \n');
        
end%switch
end
