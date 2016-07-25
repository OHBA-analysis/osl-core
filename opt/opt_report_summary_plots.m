function [opt, opt_report]=opt_report_summary_plots(opt, opt_report)

% [opt, opt_report]=opt_report_summary_plots(opt, opt_report);
%
% private fn used by opt
%
% MWW

%%%%%%%%%%%%%%%%%%%%
%% diagnostic plots over all sessions
try        
    
    if isfield(opt.results,'badchans_pre_sss') && ~all(isnan(opt.results.badchans_pre_sss))
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

    if isfield(opt.results,'bad_segments') && ~all(isnan(opt.results.bad_segments)) 
        opt_report=osl_report_set_figs(opt_report,'bad_segments');
        plot(opt.sessions_to_do,opt.results.bad_segments,'*');xlabel('sess no.');ylabel('Num bad segments');
        opt_report=osl_report_print_figs(opt_report);
    end
    
    if isfield(opt.results,'autobadoff') && ~all(isnan(opt.results.autobadoff))
        opt_report=osl_report_set_figs(opt_report,'sss_autobad_off');
        plot(opt.sessions_to_do,opt.results.autobadoff,'*');xlabel('sess no.');ylabel('SSS autobad off flag'); 
        opt_report=osl_report_print_figs(opt_report);
    end
        
    if isfield(opt.results,'icsremoved') && ~all(isnan(opt.results.icsremoved))
        opt_report=osl_report_set_figs(opt_report,'africa_ics_removed');
        plot(opt.sessions_to_do,opt.results.icsremoved,'*');xlabel('sess no.');ylabel('Num ICs removed by AFRICA'); 
        opt_report=osl_report_print_figs(opt_report);
    end
    

catch ME
    disp('Plots of summary diagnostics over sessions has failed');
    ME.getReport
end
