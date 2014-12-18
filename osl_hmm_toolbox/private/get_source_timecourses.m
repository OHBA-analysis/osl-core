function results = AB_get_source_timecourses( source_recon_results, varargin )
% Reconstruct source time courses for vector or scalar OSL beamformer
%
% varargin: strings containing additional flags:
%
% 'noreject' - do not reject bad epochs
% 'norecon'  - do not do time course reconstruction
% 'nonorm'   - do not noise normalise
% AB 2013


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
  
if any(ismember(varargin,'noreject'))
  source_recon_time_indices=1:length(source_recon_results.samples2use);
  source_recon_time_indices=source_recon_time_indices(samples_in_woi);
else
  source_recon_time_indices=find(source_recon_results.samples2use);
end


% setup channels
modality = 'MEG';
chanindmeg = strmatch(modality, D.chantype);
chanind = setdiff(chanindmeg, D.badchannels);
if isempty(chanind)
    error(['No good ' modality ' channels were found.']);
end


% load in sensor data
results.sensor_data = D(chanind, source_recon_time_indices,triallist);
results.chanind   = chanind;
results.time_inds = source_recon_time_indices;
results.triallist = triallist;
results.tbad      = epoch2logical(get_bad_sections(D),D.time(samples_in_woi));
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
  classchanind=find(strcmp(upper(D.chanlabels),'CLASS'));
  if(isempty(classchanind))
    error(['No ''CLASS'' chanlabel in: ' D.fname]);
  else
    for kk=1:NK
      class_samples_inds{kk} = (D(classchanind,source_recon_time_indices,triallist)==kk);
    end
  end
else
  class_samples_inds{1}= ones(1,length(source_recon_time_indices),length(triallist));
end

results.class_samples_inds = class_samples_inds;


%% Source recon

Nvoxels = numel(source_recon_results.BF.inverse.W.MEG{1});

% Get weights for all voxels & classes:
weights    = cell(Nvoxels,NK);
wnorm      = cell(Nvoxels,NK);
wnorm_nai  = cell(Nvoxels,NK);
variance   = cell(Nvoxels,NK);

for vox = 1:Nvoxels
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
for vox = 1:Nvoxels
  for kk = 1:NK
    wnorm{vox,kk} = diag(weights{vox,kk}*source_recon_results.BF.features.C.MEG{kk}*weights{vox,kk}');
    if ~strcmp(source_recon_results.recon_method,'none') && ~any(ismember(varargin,'nonorm'))
      % this one represents the uncertainty and will be
      % applied equally to the data and regressors, is
      % therefore irrelevant if NK=1
      % this one represents a scaling of just the data
      wnorm_nai{vox,kk} = diag(weights{vox,kk}*weights{vox,kk}');      
    else
      wnorm_nai{vox,kk} = ones(size(weights{1},1),1);
    end
      variance{vox,kk} = wnorm{vox,kk}./wnorm_nai{vox,kk};

  end
end

results.wnorm_nai = wnorm_nai;
results.weights = weights;
results.mni_coord = source_recon_results.mni_coord;
results.variance = variance;

if ~any(ismember(varargin,'norecon'))
  results.source_timecourses = get_source_timecourses_recon(results,1:Nvoxels);
end

end



function l = epoch2logical(ev,t)
l = false(size(t));

for i=1:size(ev,1)
  l = l | ( t >= ev(i,1) & t < ev(i,2) );
end

end




function badsections = get_bad_sections(D,output_format)
% output format: 'logical' or 'indexed'
% AB 2013
%
% Changed by Giles Colclough 25 Nov 2013 to account for multiple trials. If
% there is more than one trial, D.events is a cell array of structures, one
% for each trial. 

if nargin<2
  output_format = 'indexed';
end

if D.ntrials > 1,
    badsections = [];
    Events = events(D, ':');
    for iTrial = 1:D.ntrials,
        badsections = [badsections; ...
                       get_bad_sections_from_events(Events{iTrial}, ...
                                                    output_format)];       %#ok<AGROW>
    end%if
    
elseif D.ntrials == 1,
    badsections = get_bad_sections_from_events(events(D,1), output_format);
    
else
    error('What''s up with the number of trials?\n');
end%if




function badsections = get_bad_sections_from_events(Events, output_format)
if ~isempty(Events)
  Events = Events(strcmp({Events.type},'BadEpoch'));
  
  switch output_format
    case 'logical'
      badsections = false(1,D.nsamples);
      for ev = 1:numel(Events)
        badsections = badsections | D.time >= Events(ev).time & D.time < (Events(ev).time+Events(ev).duration);
      end
    case 'indexed'
        duration = [Events.duration]; 
        time     = [Events.time];
        if isempty(duration), duration = zeros(size(time)); end % catch case of duration being empty in each structure (it happens in some of the faces_subject_1 test data!)
        
        badsections = [time; time + duration]';
      
  end%switch

else
  badsections = [];
  
end%if

end

end

