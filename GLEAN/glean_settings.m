function settings = glean_settings(settings)
% !!! THIS FUNCTION IS UNDER CONSTRUCTION !!!

if isempty(settings)
    settings = struct;
end

% --- VALIDATE TOP LEVEL FIELDS --- %
p = inputParser;
p.PartialMatching = 0;
p.KeepUnmatched   = 1;
addParameter(p,'envelope', ...
               struct,    ...          
               @isstruct)
           
addParameter(p,'subspace', ...  
               struct('pca',struct), ...
               @(x) isstruct(x) ...
                 && isscalar(fieldnames(x)) ...
                 && any(strcmp(fieldnames(x),{'pca','parcellation','voxel'})))
            
addParameter(p,'model', ...
               struct('hmm',struct), ...
               @(x) isstruct(x) ...
                 && isscalar(fieldnames(x)) ...
                 && any(strcmp(fieldnames(x),{'hmm','ica'})))
              
addParameter(p,'output', ...
               struct, ...
               @isstruct)
            
parse(p,settings)
checkfields(p,settings)
settings = p.Results;
settings = orderfields(settings,{'envelope','subspace','model','output'});



% --- VALIDATE ENVELOPE FIELDS --- %
p = inputParser;
p.PartialMatching = 0;
p.KeepUnmatched   = 1;
addOptional(p,'log', ...
              0, ...
              @(x) x==0 || x == 1)
           
addOptional(p,'winsize', ...
              0.1, ...
              @(x) isnumeric(x) && isscalar(x) && (x > 0))
           
addOptional(p,'normalisation', ...
              'voxel', ...
              @ischar)
          
addOptional(p,'weights_normalisation', ...
              0, ...
              @(x) x==0 || x == 1)
          
addOptional(p,'parcellation', ...
              'none', ...
              @(x) isstruct(x) || strcmp(x,'none'))
              
parse(p,settings.envelope)
checkfields(p,settings.envelope)
settings.envelope = p.Results;

% Validate envelope.parcellation fields:
if ~strcmp(settings.envelope.parcellation,'none')
    p = inputParser;
    p.PartialMatching = 0;
    p.KeepUnmatched   = 1;
    
    addRequired(p,'file', ...
                  @ischar)
              
    addRequired(p,'mask', ...
                  @ischar)
              
    addOptional(p,'method', ...
                  'spatialBasis', ...
                  @ischar)
              
    addOptional(p,'orthogonalisation', ...
                  'none', ...
                  @ischar)
                  
    parse(p,settings.envelope.parcellation)
    checkfields(p,settings.envelope.parcellation)
    settings.envelope.parcellation = p.Results;
end



% --- VALIDATE SUBSPACE FIELDS --- %
p = inputParser;
p.PartialMatching = 0;
p.KeepUnmatched   = 1;
subspace = char(fieldnames(settings.subspace));
switch subspace
    case 'pca'
        addOptional(p,'dimensionality', ...
                      40, ...
                      @(x) isnumeric(x) && isscalar(x) && (x > 0))
                  
        addOptional(p,'whiten', ...
                      0, ...
                      @(x) x==0 || x == 1)
        
    case 'parcellation'
        addRequired(p,'file')
                  
        addOptional(p,'mask',[],@ischar) % should be addRequired, but there is a Matlab bug that breaks this :-(
                  
        addOptional(p,'method', ...
                      'spatialBasis', ...
                      @(x) any(strcmp(x,{'spatialBasis'})))
                  
    case 'voxel'
        settings.subspace.voxel = [];
        
end
parse(p,settings.subspace.(subspace))
checkfields(p,settings.subspace.(subspace))
settings.subspace.(subspace) = p.Results;

  

% --- VALIDATE MODEL FIELDS --- %
p = inputParser;
p.PartialMatching = 0;
p.KeepUnmatched   = 1;
model = char(fieldnames(settings.model));
switch char(fieldnames(settings.model))
    case 'hmm'
        addOptional(p,'nstates', ...
                      8, ...
                      @(x) isnumeric(x) && isscalar(x) && (x > 0))
                  
        addOptional(p,'nreps', ...
                      1, ...
                      @(x) isnumeric(x) && isscalar(x) && (x > 0))

    case 'ica'
        addOptional(p,'whiten', ...
                      'pca', ...
                      @(x) any(strcmp(x,{'pca','parcellation','voxel'})))
      
end
parse(p,settings.model.(model))
checkfields(p,settings.model.(model))
settings.model.(model) = p.Results;



% --- VALIDATE OUTPUT FIELDS --- %


for output = fieldnames(settings.output)'
    switch char(output)
        case {'pcorr','connectivity_profile'}
            validateattributes(settings.output.(char(output)).format,{'char'},{},mfilename,'format')
            validateattributes(settings.output.(char(output)).mask,{'char'},{},mfilename,'mask')

    end
end


end


function checkfields(p,callerID)

    extrafields = fieldnames(p.Unmatched);
    if ~isempty(extrafields)
        errorstr = 'Unrecognised field(s): ';
        for f = 1:length(extrafields)
            errorstr = [errorstr extrafields{f}];
            if f < length(extrafields)
                errorstr = [errorstr ', '];
            end
        end
        errorstr = [errorstr ' in ' inputname(2)];
        error(errorstr)
    end
    
end
