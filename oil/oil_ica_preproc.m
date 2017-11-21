function [niftiNorm, niftiNoNorm] = oil_ica_preproc(oil, reconResultsFname)
%osl_ica_preproc prepares data for ICA by enveloping and downsampling
%
%  [normResult, noNormResult] = osl_ica_preproc(OIL,RECON_RESULTS_FNAME)
%
%  [normResult, noNormResult] = osl_ica_preproc(OIL, D)
%  
%  Input contains:
%     REQUIRED: OIL.SOURCE_RECON,   containing beamforming results
%               OIL.ENVELOPING,     containing parameters for this function
%               RECON_RESULTS_FNAME source recon results file name for this
%                                   session, from a successful oat.
%          or,  D                   a source-reconstructed SPM MEEG object       
%     
%     OIL.ENVELOPING contains:
%               window_length            : the window over which the envelope is downsampled (default 1s).   
%               timewindow               : subset of trial window to remove edge effects (default is 'all').
%               ss                       : spatial smoothing of down-sampled envelopes by gaussian kernel (FWHM) (default is 4mm).               
%               gridstep                 : spatial down-sampling of data (default 8mm).
%
%  Output contains:      
%               normResult           : name of the nifti file output 
%               noNormResult         : name of the nifti file output for
%                                      the un-normalised data.

%
%  HL 060213
%  GC 15-05-2015

%% Initialise Variables

if isfield(oil.enveloping,'window_length'), window_length =  oil.enveloping.window_length; else warning('window_length not specified, Set to 1s by default');   window_length = 1;     end;
if isfield(oil.enveloping,'ss'),            ss            =  oil.enveloping.ss;            else warning('ss not specified, Set to 4mm by default');             ss            = 4;     end;
if isfield(oil.enveloping,'gridstep'),      ds            =  oil.enveloping.gridstep;      else warning('ds not specified, Set to 8mm by default');             ds            = 8;     end;
if isfield(oil.enveloping,'timewindow'),    timewindow    =  oil.enveloping.timewindow;    else warning('timewindow not specified, Using whole trial window,'); timewindow    = 'all'; end;

saveDir = fullfile(oil.source_recon.dirname, oil.enveloping.name);
ROInets.make_directory(saveDir);

%% Load Beamformed Data
try % passed in as meeg
    D = spm_eeg_load(reconResultsFname);
    [~, saveNameStem] = fileparts(D.fname);
catch ME 
    if strfind(ME.message, 'doesn''t contain SPM M/EEG data'),
        % passed in as results structure from an oat
        [~, saveNameStem]  = fileparts(reconResultsFname);
        sourceReconResults = oat_load_results(oil, saveNameStem);
        D                  = spm_eeg_load(sourceReconResults.BF.write.spmeeg.files{1});
    else
        % another error
        rethrow(ME);
    end%if
end%if

% we're expecting a particular arrangement of montages, otherwise the code
% will break
montageCheck = D.montage('getmontage', 2);
assert(any(strfind(montageCheck.name, 'with weights normalisation')),        ...
       [mfilename ':UnexpectedMontage'],                                     ...
       ['%s: Expected source-space data to have two montages, one without ', ...
        'and one with weights normalisation. \n'], mfilename);

%% Spit out variance maps as a sanity check
varianceDir  = fullfile(saveDir, 'variance-maps', filesep);
varianceFile = fullfile(varianceDir, [saveNameStem '_noise_corrected_variance_map']);
ROInets.make_directory(varianceDir);
nii.quicksave(osl_source_variance(D.montage('switch',2)), ... % second montage is weights-normalised
              varianceFile, oil.source_recon.gridstep, 2);

%% Enveloping 
% select time points to use
t = D.time;
if strcmp(timewindow, 'all'),
    timewindow = [t(1) t(end)];
else
    % timewindow should be a 2-component vector of times to use
    if ~isnumeric(timewindow) || ~isequal(length(timewindow), 2),
        error([mfilename ':WrongTimewindowFormat'], ...
              ['timewindow must be ''all'' ', ...
               'or a 2-component vector. \n']);
    end%if 
end%if

tSamples = ROInets.setdiff_pos_int(                                       ...
    find(t>=timewindow(1), 1, 'first'):find(t<=timewindow(2), 1, 'last'), ...
    ~goodsamples(D));

% create a cleaned temporary object
Dsensor     = D.montage('switch');
Dtmp        = clone(Dsensor, tempname,  ...
                    [nchannels(Dsensor), length(tSamples), D.ntrials], 0);
cleanDtmp   = onCleanup(@() delete(Dtmp));
Dtmp(:,:,:) = Dsensor(:,tSamples,:);

% do the enveloping
fprintf('Converting to time domain, enveloping and down-sampling\n')
Dtmp      = Dtmp.montage('switch', 2);
Dh_norm   = osl_hilbenv(struct('D', Dtmp, 'winsize', window_length, ...
                               'prefix', 'h'));
Dh_nonorm = osl_hilbenv(struct('D', Dtmp.montage('switch', 1), ...
                               'winsize', window_length,    ...
                               'prefix', 'hnn'));

cleanDh_norm   = onCleanup(@() delete(Dh_norm));
cleanDh_nonorm = onCleanup(@() delete(Dh_nonorm));
delete(cleanDtmp);

%% Convert to MNI nii space and save nii
try    
    fprintf('\n%s%1.1f%s%1.1f%s\n', ...
            'Saving Downsampled Envelopes as .nii files and applying ', ...
            ss, 'mm spatial smoothing and ', ds, 'mm spatial resampling.');
        
    niftiNorm   = fullfile(saveDir, [saveNameStem '_winavHE_delta_' ...
                                  num2str(window_length) 's']);
    niftiNoNorm = fullfile(saveDir, [saveNameStem '_NoWeightsNorm_winavHE_delta_' ...
                                  num2str(window_length) 's']);
                              
    nii.quicksave(unwrap_trials(Dh_norm),   niftiNorm,   oil.source_recon.gridstep);
    nii.quicksave(unwrap_trials(Dh_nonorm), niftiNoNorm, oil.source_recon.gridstep);

catch ME
    % there's a limit on the amount of data that a nifti will take
    fprintf('Error: %s\n Message: %s\n\nProgram flow maintained. \n', ...
            ME.identifier, ME.message);
    fprintf('\n%s\n', ...
            ['Unable to save as .nii. ', ...
             'Saving as .mat instead. ', ...
             'No spatial smoothing or spatial resampling will be applied.']);
    niftiNorm   = fullfile(saveDir, ...
                      [saveNameStem '_winavHE_delta_' ...
                       num2str(window_length) 's' '_ss' num2str(ss) ...
                       'mm' '_ds' num2str(ds) 'mm']);
    niftiNoNorm = fullfile(saveDir, ...
                      [saveNameStem '_NoWeightsNorm_winavHE_delta_' ...
                       num2str(window_length) 's' '_ss' num2str(ss) ...
                       'mm' '_ds' num2str(ds) 'mm']);
    
    ica_course        = unwrap_trials(Dh_norm);
    ica_course_nonorm = unwrap_trials(Dh_nonorm);
    
    save(niftiNorm,'ica_course');
    niftiNorm = fnamec;
    
    save(niftiNoNorm,'ica_course_nonorm');
    niftiNoNorm = fnamecnn;
    return
end%try

%% spatial smoothing and downsample
niftiNorm   = smooth_and_downsample(niftiNorm,   oil.source_recon.gridstep, ss, ds);
niftiNoNorm = smooth_and_downsample(niftiNoNorm, oil.source_recon.gridstep, ss, ds);
    
end%oil_ica_preproc

function data = unwrap_trials(D)
%UNWRAP_TRIALS concatenates all trial data in MEEG object D.
%
% DATA = UNWRAP_TRIALS(D)
if strcmpi(D.type, 'continuous'),
    data = D(:,:,:);
else
    trialInds = ROInets.setdiff_pos_int(1:D.nTrials, D.badtrials);
    nTrials   = length(trialInds);
    
    data = NaN(D.nchannels, D.nSamples * nTrials);
    for iTrial = 1:nTrials,
        fillInd = ((iTrial - 1) * D.nSamples + 1) : (iTrial * D.nSamples);
        data(:,fillInd) = D(:,:,trialInds(iTrial));
    end%for
end%if
end%unwrap_trials

function newName = smooth_and_downsample(fname, gs, ss, ds)
%SMOOTH_AND_DOWNSAMPLE a nifti file
%
% NEWNAME = SMOOTH_AND_DOWNSAMPLE(NAME, GS, SS, DS) smooths nifti file NAME
%   with spatial resolution GS using a kernel of width SS and downsamples
%   to DS.

OSLDIR = getenv('OSLDIR');
% Weights Normalised Data
fname2      = [fname '_ss' num2str(ss) 'mm'];
maskImage   = fullfile(OSLDIR, 'std_masks',      ...
                       ['MNI152_T1_' num2str(gs) ...
                        'mm_brain_mask']);
dsBrain     = fullfile(OSLDIR, 'std_masks',      ...
                       ['MNI152_T1_' num2str(ds) ...
                        'mm_brain.nii.gz']);
dsBrainMask = fullfile(OSLDIR, 'std_masks',    ...
                     ['MNI152_T1_' num2str(ds) ...
                      'mm_brain_mask']);

% spatial convolution
runcmd(['fslmaths ' fname     ' -s ' num2str(ss) ...
                  ' -mas ' maskImage ' tmp1']);
runcmd(['fslmaths ' maskImage ' -s ' num2str(ss) ...
                  ' -mas ' maskImage ' tmp2']);
runcmd(['fslmaths tmp1  -div tmp2 '  fname2]);
runcmd('rm tmp1.nii.gz tmp2.nii.gz');

fname3 = [fname2 '_ds' num2str(ds) 'mm'];

% Spatial Downsampling
runcmd(['flirt -in ' fname2 ' -applyxfm -init ' ...
                  getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' fname3 ...
                  ' -paddingsize 0.0 -interp trilinear -ref ' dsBrain]);      
% mask to reduce edge blurring
runcmd(['fslmaths ' fname3 ' -mas ' dsBrainMask ' ' fname3]);             

[~, nam, ext] = fileparts(fname3);
newName       = [nam ext];
end%smooth_and_downsample
% [EOF]