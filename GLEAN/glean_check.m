function GLEAN = glean_check(settings)
% !!! THIS FUNCTION IS UNDER CONSTRUCTION !!!
%
% GLEAN = glean_check(settings)

% Top level fields:
p = inputParser;
p.PartialMatching = 0;
p.KeepUnmatched   = 1;
addParameter(p,'envelope',  struct, @isstruct)
addParameter(p,'subspace',  struct('pca',struct),@(x) isstruct(x) && isscalar(fieldnames(x)) && any(strcmp(fieldnames(x),{'pca','parcellation','voxel'})))
addParameter(p,'model',     struct('hmm',struct),@(x) isstruct(x) && isscalar(fieldnames(x)) && any(strcmp(fieldnames(x),{'hmm','ica'})))
addParameter(p,'output',    struct, @isstruct)
parse(p,settings)
checkfields(p,settings)
settings = p.Results;

% Envelope settings
if ~isfield(settings.envelope,'log'), settings.envelope.log = 0; end
if ~isfield(settings.envelope,'log'), settings.envelope.log = 0; end
if ~isfield(settings.envelope,'log'), settings.envelope.log = 0; end
if ~isfield(settings.envelope,'log'), settings.envelope.log = 0; end








% Subspace settings
settings.subspace.pca.dimensionality  = 40;
settings.subspace.pca.whiten          = 1;
%settings.subspace.parcellation.file       = fullfile(DROPBOXDIR,'/Projects/HCP/HCP_parcellation_8mm.nii.gz');
%settings.subspace.parcellation.mask       = fullfile(DROPBOXDIR,'/Projects/HCP/HCP_MNI152_T1_8mm_brain_mask.nii.gz');
%settings.subspace.parcellation.method     = 'spatialbasis';

% Decompositon settings
settings.model.hmm.nstates = 8;
settings.model.hmm.nreps   = 3;
%settings.model.ica.order = 20;

% Output settings
settings.output{1}.method = 'pcorr';
settings.output{1}.format = 'nii';
settings.output{1}.mask   = fullfile(DROPBOXDIR,'/Projects/HCP/HCP_MNI152_T1_8mm_brain_mask.nii.gz');

end



function checkfields(p,~)

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