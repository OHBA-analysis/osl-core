function [opt, opt_report]=opt_report_summary_plots(opt, opt_report)

% [opt, opt_report]=opt_report_summary_plots(opt, opt_report);
%
% private fn used by opt
%
% MWW

%%%%%%%%%%%%%%%%%%%%
%% diagnostic plots over all sessions
try,
    
    opt.results.badchans=nan(length(opt.sessions_to_do),1);    
    opt.results.rejected=nan(length(opt.sessions_to_do),1);
    opt.results.autobadoff=nan(length(opt.sessions_to_do),1);
    opt.results.icsremoved=nan(length(opt.sessions_to_do),1);
    opt.results.pre_sss_badchans=nan(length(opt.sessions_to_do),1);
        
    for subi=1:length(opt.sessions_to_do), subnum=opt.sessions_to_do(subi);
        
        try,
            % load in opt results for this session:            
            opt_results=opt_load_results(opt, opt.results.fnames{subnum});
           
            S2=[];
            if opt.epoch.do, 
                spm_file = fullfile(opt.dirname, opt_results.spm_files_epoched_basename);
            else
                spm_file = fullfile(opt.dirname, opt_results.spm_files_basename);
            end

            S2.D = spm_file;   
            D=spm_eeg_load(S2.D);

            opt.results.badchans(subi)=length(D.badchannels);
            opt.results.rejected(subi)=length(D.badtrials);

            if isfield(opt_results,'maxfilter') && isfield(opt_results.maxfilter,'badchans_pre_sss')                    
                opt.results.badchans_pre_sss(subi)=length(opt_results.maxfilter.badchans_pre_sss);
            end;
            
            if isfield(opt_results,'maxfilter') && isfield(opt_results.maxfilter,'sss_autobad_off')        
                opt.results.autobadoff(subi)=opt_results.maxfilter.sss_autobad_off(subnum);
            end;

            if isfield(opt_results,'africa') && isfield(opt_results.africa, 'bad_components')
                opt.results.icsremoved(subi)=length(opt_results.africa.bad_components);
            end;
                        
        catch ME,
            disp(['Could not get summary diagnostics for ' opt.results.fnames{subnum}]);
            ME.getReport
        end;
        
    end;

    if isfield(opt_results,'maxfilter') && isfield(opt_results.maxfilter,'badchans_pre_sss')                    
        opt_report=osl_report_set_figs(opt_report,'num_bad_chans_pre_sss');
        plot(opt.sessions_to_do,opt.results.badchans_pre_sss,'*');xlabel('sess no.');ylabel('Num channels rejected pre SSS'); % plot the number of rejected channels
        opt_report=osl_report_print_figs(opt_report);
    end;
    
    opt_report=osl_report_set_figs(opt_report,'num_reject_trials');
    plot(opt.sessions_to_do,opt.results.rejected,'*');xlabel('sess no.');ylabel('Num trials rejected'); % plot the number of rejected trials        
    opt_report=osl_report_print_figs(opt_report);

    opt_report=osl_report_set_figs(opt_report,'num_bad_chans');
    plot(opt.sessions_to_do,opt.results.badchans,'*');xlabel('sess no.');ylabel('Num channels rejected'); % plot the number of rejected channels
    opt_report=osl_report_print_figs(opt_report);

    if isfield(opt_results,'maxfilter') && isfield(opt_results.maxfilter,'sss_autobad_off')                    
        opt_report=osl_report_set_figs(opt_report,'sss_autobad_off');
        plot(opt.sessions_to_do,opt.results.autobadoff,'*');xlabel('sess no.');ylabel('SSS autobad off flag'); % plot the number of rejected channels
        opt_report=osl_report_print_figs(opt_report);
    end;
        
    if isfield(opt_results,'africa')                        
        opt_report=osl_report_set_figs(opt_report,'africa_ics_removed');
        plot(opt.sessions_to_do,opt.results.icsremoved,'*');xlabel('sess no.');ylabel('Num ICs removed by AFRICA'); % plot the number of rejected channels
        opt_report=osl_report_print_figs(opt_report);
    end;
    
catch ME,
    disp('Plots of summary diagnostics over sessions has failed');
    ME.getReport
end;
