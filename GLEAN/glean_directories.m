function GLEAN = glean_directories(GLEAN)
% Sets up the directory structure for the GLEAN analysis.
% Default directories will be created if necessary based on the particular 
% GLEAN settings used.
%
% D = GLEAN_DIRECTORIES(GLEAN)
%
% If full directories aren't explicitly stated then the following directory 
% structure is created:
% 
% - ENVELOPES_DIR
%   - data
%   - SUBSPACE_DIR
%     - data
%     - MODEL_DIR
%       - model
%       - RESULTS_DIR
%           - results
%
% Adam Baker 2015

% Loop through each module and set up the directory and file structure
for module = {'envelope','subspace','model','results'}
    GLEAN = setup_dir(GLEAN,module);
end

function GLEAN = setup_dir(GLEAN,module)
    % Set up the list of directories for each stage
    if ~isempty(GLEAN.(char(module)).settings.dir)
        dirname = GLEAN.(char(module)).settings.dir;
        if dirname(end) == '/'
            dirname(end) = [];
        end
    else
        dirname = get_default_dirname(GLEAN,module);
    end
    [parentDir,subDir,~] = fileparts(dirname);

    if isempty(parentDir)
        parentDir = get_default_parentdir(GLEAN,module);
    end
    fulldir = fullfile(parentDir,subDir);
    if ~isdir(fulldir)
        mkdir(fulldir)
    end
    GLEAN.(char(module)).dir = fulldir;
end


function dirname = get_default_dirname(GLEAN,module)
% Create a default directory name for a module based on its settings
    switch char(module)

        case 'envelope'
            % envelope_R[fsample]_L[log]_F[f1l-f1h_f2l-f2h]
            dirname = sprintf('%s%s%d%s%d%s%s','envelope', ...
                             '_L',GLEAN.envelope.settings.log, ...
                             '_R',GLEAN.envelope.settings.fsample);
            if isfield(GLEAN.envelope.settings,'freqbands');
                dirname = [dirname '_F',fbandstr(GLEAN.envelope.settings.freqbands)];
            end

        case 'subspace'
            % [method]_$[method_setting1]_$[method_setting2]_..."
                switch char(intersect(fieldnames(GLEAN.subspace.settings),{'pca','parcellation','voxel'}))
                    case 'pca'
                        dirname = sprintf('%s%s%d%s%d%s%s','pca', ...
                                          '_D',GLEAN.subspace.settings.pca.dimensionality, ...
                                          '_W',GLEAN.subspace.settings.pca.whiten, ...
                                          '_N',GLEAN.subspace.settings.normalisation(1));
                    case 'parcellation'
                        [~,parcellation_fname,~] = fileparts(GLEAN.subspace.settings.parcellation.file);
                        parcellation_fname = strrep(parcellation_fname,'.nii','');
                        dirname = sprintf('%s%s%s%s%s%s',parcellation_fname,'_M',GLEAN.subspace.settings.parcellation.method(1), ...
                                          '_N',GLEAN.subspace.settings.normalisation(1));
                    case 'voxel'
                        dirname = sprintf('%s%s%s','voxel', ...
                                          '_N',GLEAN.subspace.settings.normalisation(1));
                end

        case 'model'
            model = char(intersect(lower(fieldnames(GLEAN.model.settings)),{'hmm','ica'}));
            % [method]_$[method_setting1]_$[method_setting2]_..."
            switch model
                case 'hmm'
                    dirname = sprintf('%s%s%d','hmm', ...
                                      '_K',GLEAN.model.settings.hmm.nstates);
                case 'ica'
                    dirname = sprintf('%s%s%d','ica', ...
                                       '_O',GLEAN.model.settings.ica.order);
            end

        case 'results'    
            dirname = 'results';
    end
end


function parentDir = get_default_parentdir(GLEAN,module)
% Set default parent directories
    switch char(module)
        case 'envelope'
            parentDir = fileparts(GLEAN.name);
        case 'subspace'
            parentDir = GLEAN.envelope.dir;
        case 'model'
            parentDir = GLEAN.subspace.dir;
        case 'results'
            parentDir = GLEAN.model.dir;
    end
end


function str = fbandstr(fbands)
% Convert a list of frequency bands into a single string
    tmpstr = regexp(num2str(cat(2,fbands{:})),'\s+','split');
    str = cell(1,length(tmpstr)*2 - 1);
    str(1:4:end) = tmpstr(1:2:end);
    str(3:4:end) = tmpstr(2:2:end);
    [str{2:4:end}] = deal('-');
    [str{cellfun(@isempty,str)}] = deal('_');
    str = [str{:}];

end



end
