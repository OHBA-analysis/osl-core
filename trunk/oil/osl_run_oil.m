function oil = osl_run_oil(oil)

% oil = osl_run_ica(oil)
%
% Runs OSL's ICA pipeline using the settings contained in "oil", which 
% needs to be setup by calling oil=osl_setup_ica(S). All output results are
% saved to disk and stored in the "oil" structure.
%
% Source reconstruction must have already been performed. 
%
% There are 5 stages which can be run:
%
% 1.) Estimation of down-sampled oscillatory envelope at each voxel
%     - oil.to_do(2) = 1.
%     - performed by osl_ica_preproc.m
%     - enveloping parameters and output stored in "oil.enveloping"
%     - beamformer result is converted to time domain, hilbert envelope
%       estimated and down-sampled. Envelopes are spatially smoothed and
%       resampled before being saved.
%       Analysis works for epoched data and continuous data.
% 2.) Concatentation and Normalisation
% 3.) ICA decomposition
% 4.) Subject Level Stats
% 5.) Group Level Stats
% 
% Henry Luckhoo (henry.luckhoo@trinity.ox.ac.uk)
% Version 1.1
% 010213

%% Initial admin
oil = osl_check_oil(oil);

oil.results = struct();
oil.results.date = date;

% set logfile up
oil.results.plotsdir = fullfile(oil.source_recon.dirname, 'plots');
logFileName = ['log-' datestr(now, 29) '.txt'];
ROInets.make_directory(oil.source_recon.dirname);
oil.results.logfile = fullfile(oil.results.plotsdir, logFileName);

% set first level diagnostic report up    
oil_report           = osl_report_setup(oil.results.plotsdir, ...
                                        'OIL report', ...
                                        oil.results.logfile);  
diary(oil.results.logfile);

% turns figure plotting on or off in SPM calls
spm_get_defaults('cmdline', ~oil.do_plots);


%% Conversion to Time Domain, Hilbert Enveloping, Down-sampling and Spatial Smoothing
if oil.to_do(1),
    
    if ~isfield(oil.source_recon, 'results_fnames'),
        error([mfilename ':NoSourceRecon'],                                  ...
              ['%s: Source recon stage has not been run. \n',                ...
               'To load in a previous oil, try running ',                    ...
               '"oil=osl_load_oil(oil.fname);" and then try this again.\n'], ...
              mfilename);
    end%if
       
    for iSubject = oil.source_recon.sessions_to_do
        [oil.enveloping.results.source_space_envelopes_results_fnames{iSubject},               ...
         oil.enveloping.results.source_space_envelopes_NoWeightsNorm_results_fnames{iSubject}] ...
             = oil_ica_preproc(oil, oil.source_recon.results_fnames{iSubject});
    end%for
    
    % save OIL to disk
    oil = osl_save_oil(oil);
end%if

%% Concatenation and Normalisation of Subjects
if(oil.to_do(2)), 
    
    % Check enveloping has been done
    try
        tmp = oil.enveloping.results.source_space_envelopes_results_fnames;
    catch
        error('Enveloping stage has has not been run. To load in a previous oil, try running "oil=osl_load_oil(oil.name);" and then try this again.');
    end;
    
    oil = oil_concat_subs(oil);

    % save OIL to disk
    oil = osl_save_oil(oil);
end;


%% ICA decomposition
if(oil.to_do(3)),

    if isfield(oil.concat_subs.results ,'concat_file')
        oil.ica.results.ica_concat_path{1}=[oil.source_recon.dirname '/' oil.enveloping.name '/' oil.concat_subs.name '/' oil.concat_subs.results.concat_file];
    end

    oil.ica.gridstep = oil.enveloping.gridstep;
    oil.ica.subj_ind = oil.concat_subs.results.subj_ind;
    
    oil.ica = oil_perform_ica(oil.ica);   
    
    oil.ica = rmfield(oil.ica, 'gridstep'); 
    oil.ica = rmfield(oil.ica, 'subj_ind'); 

    try
        switch lower(oil.ica.temp_or_spat)
            case {'spatial'};
                oil.ica.results.maps = oil_save_nii_ica_maps(oil,'spatial');
            case {'temporal'};
                oil.ica.results.maps = oil_save_nii_ica_maps(oil,'scaled_covariance');
        end
    catch
        warning('Unable to write ICA maps to nii file');
    end
    
    % save OIL to disk
    oil = osl_save_oil(oil);
end

%% First level stats
if(oil.to_do(4)), 
    try
        switch lower(oil.paradigm)
            case {'task'}
                tmp = oil.ica.results.mixing_matrix;
            case {'rest'}
                tmp = oil.enveloping.results.source_space_envelopes_NoWeightsNorm_results_fnames;
                tmp = oil.ica.results.mixing_matrix;
        end
    catch
        error('Enveloping or ICA stage has has not been run. To load in a previous oil, try running "oil=osl_load_oil(oil.name);" and then try this again.');
    end;
    
    switch lower(oil.paradigm)
        case {'task'}
            oil = oil_run_first_level_ica(oil);
        case {'rest'}
            oil = oil_single_subject_maps(oil);
        otherwise
            error('User must specify if data is "task" or "rest" in oil.paradigm')
    end
    
    % Save OIL to disk.
    oil = osl_save_oil(oil);
end

%% Group Level Stats
if(oil.to_do(5)), 
    try
        tmp = oil.ica_first_level.results;
    catch
        error('First level stats stage has has not been run. To load in a previous oil, try running "oil=osl_load_oil(oil.name);" and then try this again.');
    end;
    
    switch lower(oil.paradigm)
        case {'task'}
            oil = oil_run_group_level_ica(oil);
        case {'rest'}
            oil = oil_ica_maps_group_stats(oil);
        otherwise
            error('User must specify if data is "task" or "rest" in oil.paradigm')
    end
    
    % Save OIL to disk.
    oil = osl_save_oil(oil);
end

%% generate web report
oil.results.report=osl_report_write(oil_report);        

disp(['To view OIL report, point your browser to <a href="' oil.results.report.html_fname '">' oil.results.report.html_fname '</a>']);

diary off;

end%osl_run_oil
% [EOF]