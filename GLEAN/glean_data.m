function GLEAN = glean_data(GLEAN)
% Set up the directory structure for a new or existing GLEAN analysis
%
% GLEAN = glean_data(GLEAN)

% Check beamformed data exists and is in the right format
if ~isfield(GLEAN,'data') || ~all(cellfun(@ischar,GLEAN.data))
    error('Must specify GLEAN.data as a [sessions x 1] cell array"')
end
        
for module = {'envelope','subspace','model','results'}
    GLEAN = setup_dir(GLEAN,module);
    GLEAN = setup_files(GLEAN,module);
end


end



function GLEAN = setup_dir(GLEAN,module)
    if isfield(GLEAN.(char(module)).settings,'dir')
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
        % [method]_$[method_setting1]_$[method_setting2]_..."
        switch lower(char(fieldnames(GLEAN.model.settings)))
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

tmpstr = regexp(num2str(cat(2,fbands{:})),'\s+','split');
str = cell(1,length(tmpstr)*2 - 1);
str(1:4:end) = tmpstr(1:2:end);
str(3:4:end) = tmpstr(2:2:end);
[str{2:4:end}] = deal('-');
[str{cellfun(@isempty,str)}] = deal('_');
str = [str{:}];

end



function GLEAN = setup_files(GLEAN,module)

[~,sessionNames] = cellfun(@fileparts,GLEAN.data,'UniformOutput',0);

switch char(module)
    
    case {'envelope','subspace'}
        % Place data files under subdirectory "data"
        GLEAN.(char(module)).data = fullfile(GLEAN.(char(module)).dir,'data',strcat(sessionNames,'.mat'));
        if ~isdir(fullfile(GLEAN.(char(module)).dir,'data'))
            mkdir(fullfile(GLEAN.(char(module)).dir,'data'));
        end
        
    case 'model'

        GLEAN.model.model = fullfile(GLEAN.model.dir,'model.mat');
        
    case 'results'
        
        if ~isempty(fieldnames(GLEAN.results.settings))
            for field = fieldnames(GLEAN.results.settings)'
                
                results = char(field);

                mapSpaces = cellstr(GLEAN.results.settings.(results).space);
               
                for mapSpace = mapSpaces(:)'
                         
                    resultsDir = fullfile(GLEAN.results.dir,results);
                    if ~isdir(resultsDir)
                        mkdir(resultsDir);
                    end
                    
                    switch results
                        
                        case 'pcorr'
                            sessionMaps = fullfile(resultsDir,char(mapSpace),strcat(sessionNames,'_',results));
                            groupMaps   = fullfile(resultsDir,char(mapSpace),strcat('group_',results));
                            
                            % Duplicate maps across each frequency band:
                            if isfield(GLEAN.envelope.settings,'freqbands')
                                fstr      = cellfun(@(s) regexprep(num2str(s),'\s+','-'), GLEAN.envelope.settings.freqbands,'UniformOutput',0);
                                groupMaps = strcat(groupMaps,'_',fstr,'Hz.',GLEAN.results.settings.(results).format);
                                if ~isempty(sessionMaps)
                                    sessionMaps = cellfun(@(s) strcat(s,'_',fstr,'Hz.',GLEAN.results.settings.(results).format),sessionMaps,'UniformOutput',0);
                                end
                            else
                                if ~isempty(sessionMaps)
                                    sessionMaps = cellfun(@(s) {strcat(s,'.',GLEAN.results.settings.(results).format)},sessionMaps,'UniformOutput',0);
                                end
                                groupMaps = {strcat(groupMaps,'.',GLEAN.results.settings.(results).format)};
                            end
                            
                            GLEAN.results.(results).(char(mapSpace)).sessionmaps  = sessionMaps;
                            GLEAN.results.(results).(char(mapSpace)).groupmaps    = groupMaps;
                            
                        case 'connectivity_profile'
                            GLEAN.results.(results) = fullfile(resultsDir,strcat('group_',results,'.',GLEAN.results.settings.(results).format));
                    end
                end
            end
        end
end
end
