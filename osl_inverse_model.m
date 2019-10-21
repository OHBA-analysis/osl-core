function D = osl_inverse_model(D,mni_coords,varargin)
% OSL_INVERSE_MODEL runs MEG forward model in SPM12
%
% This function passes a limited set of parameters to SPM batch and returns
% an SPM object containing the inverse solution in an online montage.
%
% D = osl_inverse_model(D,mni_coords,S)
% 
% REQUIRED INPUTS:
%
% D             - SPM MEG object filename (or SPM MEG object)
%
% mni_coords    - [N x 3] list of MNI coordinates to beamform to
%
% 
% OPTIONAL INPUTS: 
%
% Passed in as a struct:
%
% S.modalities     - Sensor modalities to use (e.g. MEG,MEGPLANAR, or both) 
%                    (default MEGPLANAR for Neuromag, MEG for CTF)
% 
%
% S.type           - Beamformer type to use {'Scalar' or 'Vector'}
%                     (default 'Scalar')
%
% S.timespan       - Time range to use in seconds [start end]
%                     (default [0 Inf])
%
% S.pca_order      - PCA dimensionality to use for covariance matrix inversion
%                     (default to full rank)
%
% S.use_class_channel     
%                  - flag indicating whether or not to use the class channel
%                    in S.D to determine the time samples to be used for
%                     (default is 0)
%
% S.inverse_method - inverse method to use
%                     (default is to use 'beamform')
%
% S.conditions     - conditions to use
%                     (default is to use all)
%
% S.dirname        - dir to output results to
%                     (default is a temporary directory created within D.path)
%
% S.prefix         - write new SPM file by prepending this prefix
%                     (default is '')
%
% AB 2014, MWW 2014

arg = inputParser;
arg.addParameter('modalities',[],@(x) ischar(x) || iscell(x)); % Sensor modalities to use
arg.addParameter('fuse','no',@(x) any(strcmp(x,{'no','all','meg'}))); % fuse modalities listed in 'modalities' - e.g. set modalities to {'MEGMAG','MEGPLANAR'} to fuse them (does not work with MEGANY)
arg.addParameter('type','Scalar',@(x) any(strcmp(x,{'Scalar','Vector'})));
arg.addParameter('timespan',[0 Inf]); % Time range to use in seconds [start end]
arg.addParameter('pca_order',[],@(x) isnumeric(x) && isscalar(x) && x>0 && ~mod(x,1)); % PCA dimensionality to use for covariance matrix inversion (defaults to full rank)
arg.addParameter('use_class_channel',false); % flag indicating whether or not to use the class channel to determine the time samples to be used for
arg.addParameter('inverse_method','beamform',@(x) any(strcmp(x,{'beamform','beamform_bilateral'}))); %any(strcmp(x,{'beamform','beamform_bilateral','mne_eye','mne_diag_datacov','mne_adaptive'}))); % inverse method to use
arg.addParameter('conditions','all',@(x) ischar(x) || iscell(x)); % Should be a cell array of conditions
arg.addParameter('dirname','',@ischar); % dir to output results to - default/empty makes a temporary folder alongside the D object
arg.addParameter('prefix',''); % write new SPM file by prepending this prefix
arg.parse(varargin{:});
S = arg.Results;

%%%%%%%%%%%%%%%%%%%%%%%   P A R S E   I N P U T S   %%%%%%%%%%%%%%%%%%%%%%%

old_dir = pwd; % Back up the original working directory

if nargin < 3 || isempty(S) 
    S = struct;
end

if isa(D,'meeg')
    % If the fullfile is relative (e.g. if the user ran `D =
    % D.copy('./temp'))` then the paths fail later on. The solution is to
    % reload the file e.g. D = spm_eeg_load(D.fullfile) However, this will
    % only work properly if the user is the same directory as when they
    % created the MEEG object. We prompt the user to do this themselves,
    % rather than trying to automatically reload the file. If the user did
    % change directory e.g.
    %
    % D = D.copy('./temp')) 
    % cd .. 
    % D = spm_eeg_load(D.fullfile)
    %
    % likely scenario is that the file won't exist and cannot be loaded. Worst
    % case scenario is that the new folder contains an MEEG with the same
    % filename, in which case an incorrect file would be loaded. It should
    % definitely not be automatically reloaded
    if ~strcmp(D.fullfile,getfullpath(D.fullfile))
        error('MEEG fullfile (''%s'') is relative - reload the file (e.g. ''D = spm_eeg_load(D.fullfile)'')to make it absolute',D.fullfile);
    end
else
    if ischar(D) || isstring(D)
        [pathstr,filestr] = fileparts(D);
        D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
        D = spm_eeg_load(D);
    else
        error('Unrecognized input time - must be a file name or an MEEG');
    end
end

D.check;
D.save(); % Save the object to disk to ensure that the current online montage is used

% Check Modality Specification:
if isempty(S.modalities)
    if any(strcmp(unique(D.chantype),'MEGPLANAR'))
        default_modality = {'MEGPLANAR'};
    else
        default_modality = {'MEG'};
    end
    S.modalities = default_modality;
elseif ischar(S.modalities)
    S.modalities = {S.modalities};
end

for j = 1:length(S.modalities)
    assert(~isempty(D.indchantype(S.modalities{j},'GOOD')),'No good channels found for modality %s',S.modalities{j});
end

% Check PCA order
if isempty(S.pca_order)
    S.pca_order = D.nchannels;
    fprintf('Setting PCA order to %d (full rank)\n',D.nchannels);
end

% Check conditions Specification:
if ischar(S.conditions)
    S.conditions = {S.conditions};
end
assert(any(strcmp(S.conditions,'all')) || isempty(setdiff(S.conditions,D.condlist)),'Not all requested conditions are present in the MEEG object');

% Check dirname Specification:
if isempty(S.dirname)
    [~,tn] = fileparts(tempname(D.path));
    S.dirname = fullfile(D.path,sprintf('osl_bf_temp_%s',tn(3:3+7))); 
else
    S.dirname = getfullpath(S.dirname);
end

if ~exist(S.dirname,'dir')
    mkdir(S.dirname);
end

fprintf(1,'BF working directory: %s\n',S.dirname);

%%%%%%%%%%%%%%%%%%   R U N   I N V E R S E   M O D E L   %%%%%%%%%%%%%%%%%%

clear matlabbatch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{1}.spm.tools.beamforming.data.dir                                   = {S.dirname};
matlabbatch{1}.spm.tools.beamforming.data.D                                     = {fullfile(D.path,D.fname)};
matlabbatch{1}.spm.tools.beamforming.data.val                                   = 1;
matlabbatch{1}.spm.tools.beamforming.data.gradsource                            = 'inv';
matlabbatch{1}.spm.tools.beamforming.data.space                                 = 'MNI-aligned';
matlabbatch{1}.spm.tools.beamforming.data.overwrite                             = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{2}.spm.tools.beamforming.sources.BF(1)                              = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank                        = [2 3];
matlabbatch{2}.spm.tools.beamforming.sources.keep3d                             = 1;
matlabbatch{2}.spm.tools.beamforming.sources.visualise                          = 0;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.mni_coords.pos              = double(mni_coords);


% MESH STUFF!
%matlabbatch{2}.spm.tools.beamforming.sources.BF(1)                              = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
%matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank                        = [2 3];
%matlabbatch{2}.spm.tools.beamforming.sources.keep3d                             = 1;
%matlabbatch{2}.spm.tools.beamforming.sources.visualise                          = 0;
%matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.orient                 = 'Original';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{3}.spm.tools.beamforming.features.BF(1)                             = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));

if strcmp(S.conditions{1},'all')
    matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all            = 1;
else
    matlabbatch{3}.spm.tools.beamforming.features.whatconditions.condlabel      = S.conditions;
end

if ~S.use_class_channel
    if D.ntrials == 1
        matlabbatch{3}.spm.tools.beamforming.features.plugin.contcov            = struct([]);
    else
        matlabbatch{3}.spm.tools.beamforming.features.plugin.cov                = struct([]);
    end
else
    matlabbatch{3}.spm.tools.beamforming.features.plugin.cov_bysamples          = struct([]);
end
matlabbatch{3}.spm.tools.beamforming.features.woi                               = S.timespan*1000; % needs to be in msecs for bf_features

for jj=1:length(S.modalities),
    matlabbatch{3}.spm.tools.beamforming.features.modality{jj}                  = S.modalities{jj};
end

matlabbatch{3}.spm.tools.beamforming.features.fuse                              = S.fuse;
matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda      = 0;
matlabbatch{3}.spm.tools.beamforming.features.bootstrap                         = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabbatch{4}.spm.tools.beamforming.inverse.BF(1)                              = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));

switch S.inverse_method,
    case 'beamform'
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.pca_order     = S.pca_order;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.type          = S.type;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.bilateral     = 0;
    case 'beamform_bilateral'
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.pca_order     = S.pca_order;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.type          = S.type;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.bilateral     = 1;
    case 'mne_diag_datacov'
        mne_lambda=1;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.minimumnorm.noise_cov_type  = 'diag_datacov';
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.minimumnorm.lambda          = mne_lambda;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.minimumnorm.type            = S.type;
    case 'mne_eye'
        mne_lambda=1;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.minimumnorm.noise_cov_type  = 'eye';    
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.minimumnorm.lambda          = mne_lambda;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.minimumnorm.type            = S.type;
    case 'mne_adaptive'
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.mne_adaptive.Noise          = S.MNE.Noise;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.mne_adaptive.Options        = S.MNE.Options;
    otherwise 
        disp('Inversion method unknown!');        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{5}.spm.tools.beamforming.output.BF(1)                               = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.output.plugin.montage_osl.normalise        = 'both';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch{6}.spm.tools.beamforming.write.BF(1)                                = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg_osl.prefix             = S.prefix;

obj = onCleanup(@() cd(old_dir)); % Restore the working directory
spm_jobman('run',matlabbatch)

if ~isempty(S.prefix)
    D = spm_eeg_load(fullfile(D.path,[S.prefix D.fname]))
else
    D = spm_eeg_load(D.fullfile); % Load the file back from disk with the new online montages
end







