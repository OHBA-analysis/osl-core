%% osl_ica_preproc.m
%
%  S_out=osl_ica_preproc(OIL,RECON_RESULTS_FNAME)
%  
%  Input contains:
%     REQUIRED: OIL.SOURCE_RECON, containing beamforming results
%               OIL.ENVELOPING,   containing parameters for this function
%               RECON_RESULTS_FNAME source recon results file name for this session
%     
%     OIL.ENVELOPING contains:
%               window_length            : the window over which the envelope is downsampled (default 1s).   
%               timewindow               : subset of trial window to remove edge effects (default is 'all').
%               ss                       : spatial smoothing of down-sampled envelopes by gaussian kernel (FWHM) (default is 4mm).               
%               ds                       : spatial down-sampling of data (default 8mm).
%
%  Output structure S_out contains:      
%               fils_nifti           : name of the nifti file output 
%               fils_nifti_nonorm    : name of the nifti file output for
%                                      the un-normalised data.

%
%  HL 060213
%  Version 1.2

function S_out = oil_ica_preproc(oil, reconResultsFname)

%% Initialise Variables

if isfield(oil.enveloping,'window_length'), window_length =  S_in.window_length; else warning('window_length not specified, Set to 1s by default');   window_length = 1;     end;
if isfield(oil.enveloping,'ss'),            ss            =  S_in.ss;            else warning('ss not specified, Set to 4mm by default');             ss            = 4;     end;
if isfield(oil.enveloping,'ds'),            ds            =  S_in.ds;            else warning('ds not specified, Set to 8mm by default');             ds            = 8;     end;
if isfield(oil.enveloping,'timewindow'),    timewindow    =  S_in.timewindow;    else warning('timewindow not specified, Using whole trial window,'); timewindow    = 'all'; end;

saveDir = fullfile(oil.source_recon.dirname, oil.enveloping.name);
ROInets.make_directory(saveDir);

%% Load Beamformed Data
sourceReconResults = oat_load_results(oil,reconResultsFname);
D = sourceReconResults.BF.data.D;

%% Spit out variance maps as a sanity check
varianceDir  = fullfile(saveDir, 'variance-maps', filesep);
varianceFile = fullfile(varianceDir, [reconResultsFname '_noise_corrected_variance_map']);
ROInets.make_directory(varianceDir);
nii_quicksave(osl_source_variance(D.montage('switch',2)), ... % second montage is weights-normalised
              varianceFile, sourceReconResults.gridstep, 2);

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
    all(badsamples(D,':',':',':')));

% create a cleaned temporary object
Dtmp = copy(D, tempname);
Dtmp(:,:,:) = D(:,tSamples,:);

% do the enveloping
fprintf('Converting to time domain, enveloping and down-sampling\n')
Dtmp = Dtmp.montage('switch', 2);
Dh_norm   = osl_hilbenv(struct('D', Dtmp, 'winsize', window_length, ...
                               'prefix', 'h'));
Dh_nonorm = osl_hilbenv(struct('D', Dtmp.montage('switch', 1), ...
                               'winsize', window_length,    ...
                               'prefix', 'hnn'));
Dtmp.delete;

%% Convert to MNI nii space and save nii
try    
    fprintf('\n%s%1.1f%s%1.1f%s\n', ...
            'Saving Downsampled Envelopes as .nii files and applying ', ...
            ss, 'mm spatial smoothing and ', ds, 'mm spatial resampling.');
        
    fnamec   = fullfile(saveDir, [reconResultsFname '_winavHE_delta_' ...
                                  num2str(window_length) 's']);
    fnamecnn = fullfile(saveDir, [reconResultsFname '_NoWeightsNorm_winavHE_delta_' ...
                                  num2str(window_length) 's']);
                              
    nii_quicksave(unwrap_trials(Dh_norm),   fnamec,   sourceReconResults.gridstep);
    nii_quicksave(unwrap_trials(Dh_nonorm), fnamecnn, sourceReconResults.gridstep);

catch
    % there's a limit on the amount of data that a nifti will take
    fprintf('\n%s\n', ...
            ['Unable to save as .nii. ', ...
             'Saving as .mat instead. ', ...
             'No spatial smoothing or spatial resampling will be applied.']);
    fnamec = fullfile(saveDir, ...
                      [reconResultsFname '_winavHE_delta_' ...
                       num2str(window_length) 's' '_ss' num2str(ss) ...
                       'mm' '_ds' num2str(ds) 'mm']);
    fnamec = fullfile(saveDir, ...
                      [reconResultsFname '_NoWeightsNorm_winavHE_delta_' ...
                       num2str(window_length) 's' '_ss' num2str(ss) ...
                       'mm' '_ds' num2str(ds) 'mm']);
    
    ica_course        = unwrap_trials(Dh_norm);
    ica_course_nonorm = unwrap_trials(Dh_nonorm);
    
    save(fnamec,'ica_course');
    S_out.fils_nifti = fnamec;
    
    save(fnamecnn,'ica_course_nonorm');
    S_out.fils_nifti_nonorm = fnamec;
    return
end

%% spatial smoothing and downsample
Sout.fils.nifti        = smooth_and_downsample(fnamec,   sourceReconResults.gridstep, ss, ds);
Sout.fils.nifti_nonorm = smooth_and_downsample(fnamecnn, sourceReconResults.gridstep, ss, ds);
    
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

global OSLDIR
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
call_fsl_wrapper(['fslmaths ' fname     ' -s ' num2str(ss) ...
                  ' -mas ' maskImage ' tmp1']);
call_fsl_wrapper(['fslmaths ' maskImage ' -s ' num2str(ss) ...
                  ' -mas ' maskImage ' tmp2']);
call_fsl_wrapper(['fslmaths tmp1  -div tmp2 '  fname2]);
call_fsl_wrapper('rm tmp1.nii.gz tmp2.nii.gz');

fname3 = [fname2 '_ds' num2str(ds) 'mm'];

% Spatial Downsampling
call_fsl_wrapper(['flirt -in ' fname2 ' -applyxfm -init ' ...
                  getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' fname3 ...
                  ' -paddingsize 0.0 -interp trilinear -ref ' dsBrain]);      
% mask to reduce edge blurring
call_fsl_wrapper(['fslmaths ' fname3 ' -mas ' dsBrainMask fname3]);             

[~, nam, ext] = fileparts(fname3);
newName       = [nam ext];
end%smooth_and_downsample
% [EOF]