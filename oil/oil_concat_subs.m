%% osl_concat_subs.m
%
% osl_concat_subs is a function to be run as part of "osl_run_ica". It takes
% the down-sampled oscillatory envelope files produced by "osl_ica_preproc"
% and concatenated them into a single matrix. Optional inputs can normalise
% individual subjects (oil.concat_subs.normalise_subjects = 1) and only
% concatenate a subset of the available subjects (defined by
% oil.concat_subs.sessions_to_do). To ammend these parameter the user should edit
% the "oil" structure prior to running osl_run_ica.
%
% The output file is stored in [oil.name '/ica_dir/'].
%
% HL 060213 
% Version 1.1

function [ oil ] = oil_concat_subs(oil)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up parameters

if isfield(oil.enveloping.results,'source_space_envelopes_results_fnames');
    filnams=oil.enveloping.results.source_space_envelopes_results_fnames;
else
    error('Enveloping stage has not been run');
end

if isfield(oil.concat_subs, 'sessions_to_do');
    filnams=filnams(oil.concat_subs.sessions_to_do);
else
end

if isfield(oil.concat_subs,'normalise_subjects')
    do_gsn=oil.concat_subs.normalise_subjects;
else
    do_gsn=0;
end

if isfield(oil.concat_subs,'demean_vox')
    do_demean=oil.concat_subs.demean_vox;
else
    do_demean=1;
end

[~, filstem]=fileparts(oil.source_recon.results_fnames{oil.concat_subs.sessions_to_do(1)});
fil_out=[num2str(oil.source_recon.freq_range(1)) '-' num2str(oil.source_recon.freq_range(2))  'Hz_cat_winavHE_delta_' num2str(oil.enveloping.window_length) 's_' filstem];
if do_gsn
    fil_out=[fil_out '_gsn'];
end

if isfield(oil.enveloping,'gridstep'),     gridstep= oil.enveloping.gridstep; else     error('no spatial resolution specified'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatentate (and normalise) data

subj_glob_std=zeros(1, length(filnams));
subj_ind=[1 zeros(1, length(filnams))];
for I=1:length(filnams),           %loop to concatenate subjects
    
    try
        xmat=load([oil.source_recon.dirname '/' oil.enveloping.name '/' filnams{I}],'-mat');
        xmat=xmat.ica_course;
    catch
        xmat = nii.quickread([oil.source_recon.dirname '/' oil.enveloping.name '/' filnams{I}],gridstep);
    end
    
    subj_ind(I+1)=subj_ind(I)+size(xmat,2);
    
    if(do_gsn)
        subj_glob_std(I)=std(squash(demean(xmat,2)));
    else
        subj_glob_std(I)=1;
    end
    
    if(do_demean)
        subj_mean=repmat(mean(xmat,2),1,size(xmat,2));
    else
        subj_mean=zeros(size(xmat));
    end
    
    if(I>1)
        xall=cat(2, xall, (xmat-subj_mean)/subj_glob_std(I));
    else
        xall=(xmat-subj_mean)/subj_glob_std(I);
    end;
end;

ica_concat=xall;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save concatenated data
concatSaveDir = fullfile(oil.source_recon.dirname, oil.enveloping.name, ...
                         oil.concat_subs.name, filesep);
[~, concatFileName] = fileparts(fil_out); % cut off extension

if ~isdir(concatSaveDir),
    mkdir(concatSaveDir); 
end%if

try
    if size(ica_concat,2) < 32767
        nii.quicksave(ica_concat, ...
                      fullfile(concatSaveDir, fil_out), ...
                      gridstep);
    else
        disp('Nii file limit exceeded, saving as .mat \n')
        fil_out = [concatFileName '.mat']; % overwrite previous extension
        save(fullfile(concatSaveDir, fil_out), 'ica_concat','-v7.3');
    end%if
catch ME
    disp('Not able to save as Nii, saving as .mat \n')
    warning('Trying to save as Nii gave error: \n %s', ME.message);
    fil_out = [concatFileName '.mat']; % overwrite previous extension
    save(fullfile(concatSaveDir, fil_out), 'ica_concat','-v7.3');
end%try

oil.concat_subs.results.concat_file = fil_out;
oil.concat_subs.results.subj_ind    = subj_ind;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatentate non-weights-normalised data and save.

if isfield(oil.enveloping.results,'source_space_envelopes_NoWeightsNorm_results_fnames');
    filnams=oil.enveloping.results.source_space_envelopes_NoWeightsNorm_results_fnames(oil.concat_subs.sessions_to_do);

    for I=1:length(filnams),           %loop to concatenate subjects
        try
            xmat=load([oil.source_recon.dirname '/' oil.enveloping.name '/' filnams{I}],'-mat');
            xmat=xmat.ica_course_nonorm;
        catch
            xmat = nii.quickread([oil.source_recon.dirname '/' oil.enveloping.name '/' filnams{I}],gridstep);
        end
        
        if(I>1)
            xall=cat(2, xall, xmat);
        else
            xall=xmat;
        end;
        
        
    end;
    
    ica_concat_nonorm=xall;
    
    % save 
    [~, concatFileNameNoNorm] = fileparts(fil_out); % remove extension
    concatFileNameNoNorm = [concatFileNameNoNorm '_NoWeightsNorm'];
    try
        if size(ica_concat,2) < 32767
            fil_out = concatFileNameNoNorm;
            nii.quicksave(ica_concat_nonorm, ...
                          fullfile(concatSaveDir, fil_out), ...
                          gridstep);
        else
            disp('Nii file limit exceeded, saving as .mat \n')
            fil_out = [concatFileNameNoNorm '.mat']; 
            save(fullfile(concatSaveDir, fil_out), 'ica_concat','-v7.3');
        end%if
    catch ME2
        disp('Not able to save as Nii, saving as .mat \n')
        warning('Trying to save as Nii gave error: \n %s', ME2.message);
        fil_out = [concatFileNameNoNorm '.mat']; 
        save(fullfile(concatSaveDir, fil_out), 'ica_concat','-v7.3');
    end%try
    
    oil.concat_subs.results.concat_file_nonorm = fil_out;
    
else
    warning('Non-weights normalised data not found... ICA first level stats will fail.');
end
