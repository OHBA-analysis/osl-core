% Example beamforming script for HCP motor task data in OSL2
%
% Adam Baker 2015


% Generate a list of HCP subjects and sessions
dataDir  = '/Users/abaker/Data/HCP/';
savefile = '/Users/abaker/Copy/Projects/HCP/HCPmotor.mat';
BFdir = '/Users/abaker/Data/HCPmotor/';

% Get subject identifier:
subjIDs = dir(dataDir);
subjIDs = setdiff({subjIDs.name},{'.','..','.DS_Store'});


% Get list of all sessions:
datafiles      = {};
noisefiles     = {};
trialinfo      = {};
headmodels     = {};
sourcemodels   = {};
session_ID     = [];
session_subjID = [];

% Get info for each session:
for subjID = subjIDs
    
    % Data filename:
    sessions_subj  = osl_filelist(fullfile(dataDir,subjID,'MEG/Motort/tmegpreproc/'),'*TEMG*');
    datafiles      = [datafiles; sessions_subj];
    
    % Noisedata filename:
    noise_subj = repmat(fullfile(dataDir,subjID,'unprocessed/2-Pnoise/4D/c,rfDC'),size(sessions_subj));
    noisefiles = [noisefiles; noise_subj];
    
    % Trialinfo filename:
    trlinfo   = osl_filelist(fullfile(dataDir,subjID,'MEG/Motort/tmegpreproc/'),'*trialinfo*');
    trialinfo = [trialinfo; trlinfo];
        
    % Headmodel filename:
    headmodel_subj = osl_filelist(fullfile(dataDir,subjID,'MEG/anatomy/'),'*headmodel*');
    headmodels     = [headmodels; repmat(headmodel_subj,length(sessions_subj),1)];
    
    % Sourcemodel filename:
    sourcemodel_subj = osl_filelist(fullfile(dataDir,subjID,'MEG/anatomy/'),'*sourcemodel_3d8mm*');
    sourcemodels     = [sourcemodels; repmat(sourcemodel_subj,length(sessions_subj),1)];
    
   % Subject identifier:
    session_subjID = [session_subjID; repmat(str2num(subjID{:}),length(sessions_subj),1)];
    
    % Session identifier:
    tempstr = sessions_subj;
    pathstr = fileparts(sessions_subj{1});
    tempstr = strrep(tempstr,pathstr,'');   
    tempstr = strrep(tempstr,'-Motort_tmegpreproc_TEMG.mat','');
    tempstr = strrep(tempstr,['/' subjID{:} '_MEG_'],'');
    session_ID = [session_ID; cellfun(@str2num,tempstr)];
end

session_subjID = num2cell(session_subjID');
session_ID     = num2cell(session_ID');

% Write Dataset structure with entry for each session:
Dataset = repmat(struct,numel(datafiles),1);
[Dataset.processed]     = datafiles{:};
[Dataset.noise]         = noisefiles{:};
[Dataset.trialinfo]     = trialinfo{:};
[Dataset.headmodel]     = headmodels{:};
[Dataset.sourcemodel]   = sourcemodels{:};
[Dataset.subjectID]     = session_subjID{:};
[Dataset.sessionID]     = session_ID{:};


if ~isdir(BFdir)
    mkdir(BFdir);
end

[pathstr,filestr] = cellfun(@fileparts,{Dataset.processed},'uniformoutput',0);
beamformed_files = fullfile(BFdir,strcat(filestr,'.mat'))';

[Dataset.beamformed] = deal(beamformed_files{:});

session = 1;

% Load data:
load(Dataset(session).processed);   % loads data
load(Dataset(session).headmodel);   % loads headmodel
load(Dataset(session).sourcemodel); % loads sourcemodel3d
load(Dataset(session).trialinfo);   % loads trlInfo
freqband = [13 30];

% Make subselection of trials and channels 
chan_inds = ~cellfun(@isempty,regexp(data.label,'A\d+'));
trl_inds  = data.trialinfo(:,2) < 6; % excludes rest periods 
data.label = data.label(chan_inds);
data.trial = data.trial(trl_inds);
data.time = data.time(trl_inds);
data.trialinfo = data.trialinfo(trl_inds,:);
trlInfo.lockTrl{strcmp(trlInfo.lockNames,'TEMG')}(:,2);
for trl = 1:numel(data.trial)
    data.trial{trl} = data.trial{trl}(chan_inds,:);
end

% apply 3rd order gradiometers - copied from hcp_icamne L100.
% grad = ft_apply_montage(grad, grad.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');
headmodel     = ft_convert_units(headmodel,'mm');
data.grad     = ft_convert_units(data.grad,'mm');


% Load Raw Noise data
cfg = [];
cfg.channel ='MEG';
cfg.dataset = Dataset(session).noise;
noisedata = ft_preprocessing(cfg);
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
    data.trial{trl}  = ft_preproc_bandpassfilter(data.trial{trl}, data.fsample, freqband, 5);
end
noisedata.trial{1} = ft_preproc_bandpassfilter(noisedata.trial{1}, noisedata.fsample, freqband, 5);

% Compute the beamformer weights
C = 0;
for trl = 1:numel(data.trial)
    C  = C + osl_cov(data.trial{trl});    
end
C = C ./ numel(data.trial);

Cn = osl_cov(noisedata.trial{1}(noise_chaninds,:));
%Cn = eye(size(C));



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


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out SPM MEEG object to store the results in

D = spm_eeg_ft2spm(data,Dataset(session).beamformed);

% Condition labels:
trlLbls = cell(size(data.trial));
trlLbls(data.trialinfo(:,2)==1) = deal({'left hand'});
trlLbls(data.trialinfo(:,2)==2) = deal({'left foot'});
trlLbls(data.trialinfo(:,2)==4) = deal({'right hand'});
trlLbls(data.trialinfo(:,2)==5) = deal({'right foot'});

D = conditions(D,1:D.ntrials,trlLbls);
D.save;


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

%%
V = osl_source_variance(D);
HCP_savenii(mean(V,2),Dataset(session).sourcemodel,'/Users/abaker/Data/HCPmotor/variance.nii.gz');
fslview('/Users/abaker/Data/HCPmotor/variance.nii.gz','addmask')
%%
con = {[1 0 -1 0]',[0 1 0 -1]'};

x = double(cell2mat(arrayfun(@(x) strcmp(D.conditions,D.condlist{x})',1:4,'uniformoutput',0)));
cope    = zeros(size(V,1),numel(con));
varcope = zeros(size(V,1),numel(con));
pinvxtx = pinv(x'*x);
pinvx   = pinv(x);
for voxel = 1:size(V,1)
    [cope(voxel,:), varcope(voxel,:)] = glm_fast_for_meg(V(voxel,:)',x,pinvxtx,pinvx,con,0);
end

tstat = cope ./ sqrt(varcope);

HCP_savenii(tstat,Dataset(session).sourcemodel,'/Users/abaker/Data/HCPmotor/tstat.nii.gz');
fslview('/Users/abaker/Data/HCPmotor/tstat.nii.gz','addmask')