function opt = opt_gather_results( opt )

% opt = opt_gather_results( opt )
%
% private func used by opt_consolidate results and osl_run_opt
%
% MWW 2016

num_sessions=length(opt.convert.spm_files_basenames);

opt.results.spm_files=cell(num_sessions,1);
opt.results.spm_files_basenames=cell(num_sessions,1);
opt.results.spm_files_epoched=cell(num_sessions,1);
opt.results.spm_files_epoched_basenames=cell(num_sessions,1);
opt.results.badchans=nan(length(opt.sessions_to_do),1);    
opt.results.bad_segments=nan(length(opt.sessions_to_do),1);
opt.results.rejected=nan(length(opt.sessions_to_do),1);
opt.results.autobadoff=nan(length(opt.sessions_to_do),1);
opt.results.icsremoved=nan(length(opt.sessions_to_do),1);
opt.results.pre_sss_badchans=nan(length(opt.sessions_to_do),1);
opt.results.pre_sss_badchans=nan(length(opt.sessions_to_do),1);
opt.results.ica_files=cell(length(opt.sessions_to_do),1);
opt.results.spm_files_ica_input=cell(length(opt.sessions_to_do),1);

for subi=1:length(opt.sessions_to_do), 
    subnum=opt.sessions_to_do(subi);

    try
        % load in opt results for this session:            
        opt_results=opt_load_results(opt, opt.results.fnames{subnum});

        % spm files
        opt.results.spm_files_basenames{subnum}=opt_results.spm_files_basename;
        opt.results.spm_files{subnum}=[opt.dirname '/'  opt.results.spm_files_basenames{subnum}];

        if opt.epoch.do,

            opt.results.spm_files_epoched_basenames{subnum}=opt_results.spm_files_epoched_basename;
            opt.results.spm_files_epoched{subnum}=[opt.dirname '/'  opt.results.spm_files_epoched_basenames{subnum}];

        end

        % ica results
        if opt.africa.todo.ica || opt.africa.todo.ident || opt.africa.todo.remove
            try
                ica_S=load(opt_results.ica_file);                
            catch
                try
                    % try this:
                    opt_results.ica_file=[opt.dirname '/africa/' opt.results.spm_files_basenames{subnum} '.mat']; opt_save_results(opt, opt_results);
                    ica_S=load(opt_results.ica_file);
                catch
                    warning(['Could not load ica_file for ' opt.results.fnames{subnum}]);
                end
            end
            
            try 
                opt.results.ica_files{subnum}=ica_S.S.ica_file;
                opt.results.spm_files_ica_input{subnum}=ica_S.S.D;
            catch
                warning(['Invalid ica_file for ' opt.results.fnames{subnum}]);
            end    
        end
        
        S2=[];
        if opt.epoch.do, 
            spm_file = fullfile(opt.dirname, opt_results.spm_files_epoched_basename);
        else
            spm_file = fullfile(opt.dirname, opt_results.spm_files_basename);
        end

        S2.D = spm_file;   
        D=spm_eeg_load(S2.D);

        opt.results.badchans(subnum)=length(D.badchannels);
        opt.results.rejected(subnum)=length(D.badtrials);
        
        if isfield(opt_results,'bad_segments') && isfield(opt_results.bad_segments, 'bad_segments')
            opt.results.bad_segments(subnum)=length(opt_results.bad_segments.bad_segments);
        end
        
        if isfield(opt_results,'maxfilter') && isfield(opt_results.maxfilter,'badchans_pre_sss')                    
            opt.results.badchans_pre_sss(subnum)=length(opt_results.maxfilter.badchans_pre_sss);
        end

        if isfield(opt_results,'maxfilter') && isfield(opt_results.maxfilter,'sss_autobad_off')        
            opt.results.autobadoff(subnum)=opt_results.maxfilter.sss_autobad_off(subnum);
        end

        if isfield(opt_results,'africa') && isfield(opt_results.africa, 'bad_components')
            opt.results.icsremoved(subnum)=length(opt_results.africa.bad_components);
        end

    catch ME,
        disp(['Could not get results for subject number ' num2str(subnum)]);
        ME.getReport
    end;        

end

end

