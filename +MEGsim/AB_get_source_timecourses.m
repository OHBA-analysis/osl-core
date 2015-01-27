function results = AB_get_source_timecourses(source_recon_results, varargin)
% Reconstruct source time courses for vector or scalar OSL beamformer
%
% varargin: strings containing additional flags:
%
% 'noreject' - do not reject bad epochs
% 'norecon'  - do not do time course reconstruction
% 'nonorm'   - do not noise normalise
% 'voxelind' - list of voxel indices at which to reconstruct. 
% AB 2013

inputCharItems = cellfun(@(C) ischar(C), varargin);
inputFlags     = varargin(inputCharItems);

% load data 
D=osl_get_oat_sensordata(source_recon_results);

% setup trial list
triallist = D.pickconditions(source_recon_results.source_recon.conditions);

% setuptime indices
t = source_recon_results.BF.data.D.time;
if ~isempty(source_recon_results.source_recon.time_range)
  samples_in_woi = epoch2logical(source_recon_results.source_recon.time_range,t); % fix because woi lumped with bad sections
else
  samples_in_woi = true(size(t));
end
  
if any(ismember(inputFlags,'noreject'))
  source_recon_time_indices=1:length(source_recon_results.samples2use);
  source_recon_time_indices=source_recon_time_indices(samples_in_woi);
else
  source_recon_time_indices=find(source_recon_results.samples2use(samples_in_woi));
end


% setup channels
modality   = 'MEG';
chanindmeg = MEGsim.megchannels(D, modality);
chanind    = setdiff(chanindmeg, D.badchannels);
if isempty(chanind)
    error(['No good ' modality ' channels were found.']);
end


% load in sensor data
results.sensor_data = D(chanind, source_recon_time_indices,triallist);
results.chanind   = chanind;
results.time_inds = source_recon_time_indices;
results.triallist = triallist;
results.tbad      = epoch2logical(MEGsim.AB_get_bad_sections(D),D.time(samples_in_woi));
results.time      = D.time(samples_in_woi);
results.fsample   = D.fsample;


% OSL1.5 hacks:
if ~isfield(source_recon_results.BF.inverse,'W')
  if isfield(source_recon_results.BF.inverse.MEG.class{1},'W')
    source_recon_results.BF.inverse.W.MEG{1} = source_recon_results.BF.inverse.MEG.class{1}.W;
  end
end
if ~isfield(source_recon_results.BF.features,'C')
  if isfield(source_recon_results.BF.features.MEG.class{1},'C')
    source_recon_results.BF.features.C.MEG{1} = source_recon_results.BF.features.MEG.class{1}.C;
  end
end




% HMM stuff
NK=numel(source_recon_results.BF.inverse.W.MEG);
if NK>1
  classchanind=find(strcmpi(D.chanlabels,'CLASS'));
  if(isempty(classchanind))
    error(['No ''CLASS'' chanlabel in: ' D.fname]);
  else
    for kk=NK:-1:1
      class_samples_inds{kk} = (D(classchanind,source_recon_time_indices,triallist)==kk);
    end
  end
else
  class_samples_inds{1}= ones(1,length(source_recon_time_indices),length(triallist));
end

results.class_samples_inds = class_samples_inds;


%% Source recon
% parse input voxel list
voxelListSuppliedFlagInd = ismember(inputFlags, 'voxelind');
maxNvoxels               = numel(source_recon_results.BF.inverse.W.MEG{1});

if any(voxelListSuppliedFlagInd),
    voxelList  = varargin{find(voxelListSuppliedFlagInd) + 1};
    validateattributes(voxelList,                         ...
                       {'numeric'},                       ...
                       {'vector', 'positive', 'nonempty', ...
                        '<=', maxNvoxels},                ...
                       mfilename,                         ...
                       'voxelList');
    voxelList = unique(voxelList);
else
    voxelList = (1:maxNvoxels).';
end%if

% Get weights for all voxels & classes:
weights    = cell(maxNvoxels,NK);
wnorm      = cell(maxNvoxels,NK);
wnorm_nai  = cell(maxNvoxels,NK);
variance   = cell(maxNvoxels,NK);

for vox = 1:maxNvoxels
  for kk = 1:NK
    if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis
      weights{vox,kk}=zeros(1,length(chanind));
      % for sensor space - just compute for all sensors
      weights{vox,kk}(vox)=1;
    else
      weights{vox,kk}=source_recon_results.BF.inverse.W.MEG{kk}{vox};
    end
  end
end


% Get normalisation for all voxels & classes:
for vox = 1:maxNvoxels
  for kk = 1:NK
    if ~strcmp(source_recon_results.recon_method,'none') && ~any(ismember(inputFlags,'nonorm'))
      % this one represents the uncertainty and will be
      % applied equally to the data and regressors, is
      % therefore irrelevant if NK=1
      wnorm{vox,kk} = diag(weights{vox,kk}*source_recon_results.BF.features.C.MEG{kk}*weights{vox,kk}');
      % this one represents a scaling of just the data
      wnorm_nai{vox,kk} = diag(weights{vox,kk}*weights{vox,kk}');      
    else
      wnorm{vox,kk}     = ones(size(weights{1},1),1);
      wnorm_nai{vox,kk} = ones(size(weights{1},1),1);
    end
      variance{vox,kk} = wnorm{vox,kk}./wnorm_nai{vox,kk};
  end
end



results.wnorm = wnorm;
results.wnorm_nai = wnorm_nai;
results.weights = weights;
results.mni_coord = source_recon_results.mni_coord;
results.variance = variance;

if ~any(ismember(inputFlags,'norecon'))
  results.source_timecourses = MEGsim.AB_get_source_timecourses_recon(results, voxelList);
end

end



function l = epoch2logical(ev,t)
l = false(size(t));

for i=1:size(ev,1)
  l = l | ( t >= ev(i,1) & t < ev(i,2) );
end

end