function D = HCP_beamform(S)
% D = HCP_beamform(S)
% 
% S.megdata         - Filename specifying the preprocessed MEG data
% S.headmodel       - Filename specifying the headmodel
% S.sourcemodel     - Filename specifying the sourcemodel
% S.savefile        - Filename specfiying the path and name of the created SPM object
%
% S.freqband        - [low high]
% S.combinetrials   - Concatenate trials? [0/1] 
% 
% Adam Baker 2015


if ~isfield(S,'combinetrials')
    S.combinetrials = 1;
end


% Load data:
load(S.megdata);     % loads data
load(S.headmodel);   % loads headmodel
load(S.sourcemodel); % loads sourcemodel3d

% Make subselection of channels 
chan_inds = ~cellfun(@isempty,regexp(data.label,'A\d+'));
data.label = data.label(chan_inds);
for trl = 1:numel(data.trial)
    data.trial{trl} = data.trial{trl}(chan_inds,:);
end

% apply 3rd order gradiometers - copied from hcp_icamne L100.
% grad = ft_apply_montage(grad, grad.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');
headmodel     = ft_convert_units(headmodel,'mm');
data.grad     = ft_convert_units(data.grad,'mm');


% Concatenate trials
if S.combinetrials == 1
    
    % Reshape data by concatenating trials
    % (apply 2nd order detrending first to remove jumps)
    cfg             = [];
    cfg.polyremoval = 'yes';
    cfg.polyorder   = 2;
    data            = ft_preprocessing(cfg,data);
    data.trial      = {cat(2, data.trial{:})};
    data.time       = {(1:length([data.time{:}]))./data.fsample}; % not real time, due to bad trials having been removed
    data.trialinfo  = [nan 0];
end


% Load Raw Noise data
cfg = [];
cfg.channel ='MEG';
cfg.dataset = S.noisedata;
noisedata = ft_preprocessing(cfg);

% Downsample noise data
cfg = [];
cfg.resamplefs = data.fsample;
noisedata = ft_resampledata(cfg,noisedata);

% Select noise channels
[~,noise_chaninds] = ismember(data.label,noisedata.label);

% Compute the leadfields
cfg             = [];
cfg.vol         = headmodel;
cfg.grid        = sourcemodel3d;
cfg.grad        = data.grad;
cfg.channel     = data.label;
cfg.normalize   = 'no';
cfg.reducerank  = 2;

% prepare leadfields
gridLF = ft_prepare_leadfield(cfg);


% Bandpass filter before source localization
for trl = 1:numel(data.trial)
    data.trial{trl}  = ft_preproc_bandpassfilter(data.trial{trl}, data.fsample, S.freqband, 5);
end
noisedata.trial{1} = ft_preproc_bandpassfilter(noisedata.trial{1}, noisedata.fsample, S.freqband, 5);

% Compute the beamformer weights
C = 0;
for trl = 1:numel(data.trial)
    C  = C + osl_cov(data.trial{trl});    
end
C = C ./ numel(data.trial);
Cn = osl_cov(noisedata.trial{1}(noise_chaninds,:));


% Estimate rank of data covariance matrix:
eigDiff = diff(log(svd(C)));
rankVec = eigDiff<10*median(eigDiff);
rankVec(1:200) = false;
rankC = min([find(rankVec,1,'first'),rank(C)]);

invC = pinv_plus(C,rankC); 

% rank for leadfields
rankLF = 2;


% Compute weights for each voxel
weights = cell(size(sourcemodel3d.inside));
for voxel = 1:length(weights)
    
    lf  = gridLF.leadfield{sourcemodel3d.inside(voxel)};
    
    % Scalar beamformer - reduce leadfield to rank 1
    [u, ~] = svd(real(pinv_plus(lf' * invC *lf, rankLF, 0)),'econ');
    eta = u(:,1);
    lf = lf * eta;
    
    % LCMV weights
    weights{voxel} = pinv_plus(lf' * invC * lf, rankLF, 0) * lf' * invC;
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out SPM MEEG object to store the results in

D = spm_eeg_ft2spm(data,S.savefile);
D.save;

% Downsample
S = [];
S.D             = fullfile(D.path,D.fname);
S.fsample_new   = 250;
osl_spmfun(@spm_eeg_downsample,S);

% Set up SPM montage
montage = [];
montage.labelorg = D.chanlabels;
montage.labelnew = cell(numel(weights),1);

tra{1} = zeros(numel(weights),size(weights{1}(1,:),2));
tra{2} = zeros(numel(weights),size(weights{1}(1,:),2));

montage_name = {'without weights normalisation','with weights normalisation'};

for voxel = 1:numel(weights)
    montage.labelnew{voxel} = mat2str(sourcemodel3d.pos(voxel,:),3);
    tra{1}(voxel,:)     = weights{voxel};
    tra{2}(voxel,:)     = weights{voxel}./(sqrt(weights{voxel}*Cn*weights{voxel}'));
end

montage.chantypenew = repmat({'VE'}, length(montage.labelnew), 1);
montage.chanunitnew = repmat({'nA*m'}, length(montage.labelnew), 1);
montage.chantypeorg = chantype(D,D.indchannel(montage.labelorg))';
montage.chanunitorg = units(D, D.indchannel(montage.labelorg))';

% Online montage needs additional channel information
for ch = 1:length(montage.labelnew)
    montage.channels(ch).label = montage.labelnew{ch};
    montage.channels(ch).type  = montage.chantypenew{ch};
    montage.channels(ch).units = montage.chanunitnew{ch};
    montage.channels(ch).bad   = 0;
end

S = [];
S.D            = fullfile(D.path,D.fname);
S.montage      = montage;
S.keepsensors  = false;
S.keepothers   = false;
S.mode         = 'switch';

for m = 1:2
    S.montage.name = montage_name{m};
    S.montage.tra = tra{m};
    D = spm_eeg_montage(S);
    D.save;
end

clear sourcedata

end

