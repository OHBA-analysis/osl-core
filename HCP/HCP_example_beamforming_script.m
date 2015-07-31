% Example beamforming script for HCP resting data in OSL2
%
% Adam Baker 2015


% Generate a list of HCP subjects and sessions
dataDir  = '/Volumes/Data/HCP/';
BFdir    = '/Volumes/Data/HCP_beamformed/';

% Get subject identifier:
subjIDs = dir(dataDir);
subjIDs = setdiff({subjIDs.name},{'.','..','.DS_Store'});

% Get list of all sessions:
datafiles      = {};
noisefiles     = {};
headmodels     = {};
sourcemodels   = {};
session_ID     = [];
session_subjID = [];

% Get info for each session:
for subjID = subjIDs
    
    % Data filename:
    sessions_subj  = osl_filelist(fullfile(dataDir,subjID,'MEG/Restin/rmegpreproc/'),'*MEG*');
    datafiles = [datafiles; sessions_subj];
    
    % Noisedata filename:
    noise_subj = repmat(fullfile(dataDir,subjID,'unprocessed/MEG/2-Pnoise/4D/c,rfDC'),size(sessions_subj));
    noisefiles = [noisefiles; noise_subj];
    
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
    tempstr = strrep(tempstr,'-Restin_rmegpreproc.mat','');
    tempstr = strrep(tempstr,['/' subjID{:} '_MEG_'],'');
    session_ID = [session_ID; cellfun(@str2num,tempstr)];
end

session_subjID = num2cell(session_subjID');
session_ID     = num2cell(session_ID');

% Write Dataset structure with entry for each session:
Dataset = repmat(struct,numel(datafiles),1);
[Dataset.processed]    = datafiles{:};
[Dataset.noise]        = noisefiles{:};
[Dataset.headmodel]    = headmodels{:};
[Dataset.sourcemodel]  = sourcemodels{:};
[Dataset.subjectID]    = session_subjID{:};
[Dataset.sessionID]    = session_ID{:};

[pathstr,filestr] = cellfun(@fileparts,{Dataset.processed},'uniformoutput',0);
beamformed_files = fullfile(BFdir,strcat(filestr,'.mat'))';

[Dataset.beamformed] = deal(beamformed_files{:});


%% Run the beamformer:

% Select only first session for now:
Dataset = Dataset([Dataset.sessionID] == 3);

for session = 1:numel(Dataset)
    
    S.megdata     = Dataset(session).processed; 
    S.noisedata   = Dataset(session).noise;
    S.headmodel   = Dataset(session).headmodel; 
    S.sourcemodel = Dataset(session).sourcemodel; 
    S.savefile    = Dataset(session).beamformed; 
    S.freqband    = [4 30];
   
    HCP_beamform(S)
end
