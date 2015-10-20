function GLEAN = glean_data(GLEAN)
% Set up the directory structure for a new or existing GLEAN analysis.
%
% GLEAN = glean_data(GLEAN)
%
% Adam Baker 2015


% Check data exists and is in the right format
if ~isfield(GLEAN,'data') || ~all(cellfun(@ischar,GLEAN.data))
    error('Must specify GLEAN.data as a [sessions x 1] cell array"')
end
   
% Loop through each module and set up the directory and file structure
for module = {'envelope','subspace','model','results'}
    GLEAN = setup_files(GLEAN,module);
end

end


function GLEAN = setup_files(GLEAN,module)
% Set up the list of files for each stage
    [~,sessionNames] = cellfun(@fileparts,GLEAN.data,'UniformOutput',0);

    switch char(module)

        case {'envelope','subspace'}
            GLEAN.(char(module)).data = fullfile(GLEAN.(char(module)).dir,'data',strcat(sessionNames,'.mat'));
            if ~isdir(fullfile(GLEAN.(char(module)).dir,'data'))
                mkdir(fullfile(GLEAN.(char(module)).dir,'data'));
            end

        case 'model'
            GLEAN.model.model = fullfile(GLEAN.model.dir,'model.mat');

        case 'results'
            result_types = setdiff(fieldnames(GLEAN.results.settings),'dir')';
            if ~isempty(result_types)
                for result_type = result_types
                    
                    results = char(result_type);
                    
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




