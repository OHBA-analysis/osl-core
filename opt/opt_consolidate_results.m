function [ opt ] = opt_consolidate_results( opt )

% [ opt ] = opt_consolidate_results( opt )
%
% Searches the opt directory for files with the appropriate names
% (found in opt.convert.spm_files_basenames)
% and fills up what it can of opt.results for those files 
% that it can find.
%
% Note that this overwrites any existing opt.results
%
% MWW 2014

opt=osl_check_opt(opt);

%%%%%%%%%%%%%%%%%%%%
%% setup
% set logfile up
opt.results=[];
opt.results.plotsdir=[opt.dirname '/plots'];
opt.results.logfile=[opt.results.plotsdir '/log-' date '.txt'];
% delete any existing diary file with the same name
runcmd(['rm -f ' opt.results.logfile]);

opt.results.date=date;

% set diagnostic report up    
opt_report=osl_report_setup(opt.results.plotsdir,['OPT report'],opt.results.logfile);  
diary(opt.results.logfile);

num_sessions=length(opt.convert.spm_files_basenames);

%%%%%%%%%%%%%%%%%%%%
%% find SPM MEEG file results

for subi=1:length(opt.sessions_to_do), 
    subnum=opt.sessions_to_do(subi);
         
    % try and load in results
    opt.results.fnames{subnum}=['session' num2str(subnum)];
    
    try
        opt_results=opt_load_results(opt, opt.results.fnames{subnum});
    catch
        opt.results.fnames{subnum}=[];
    end
    
    % build logfile
    runcmd(['cat ' opt_results.logfile '>' opt.results.logfile]);
        
end

% gather the results together
opt = opt_gather_results(opt);

%%%%%%%%%%%%%%%%%%%%
%% do web report

% write links to sub reports written previously
for subi=1:length(opt.sessions_to_do), subnum=opt.sessions_to_do(subi);
    opt_report.sub_reports{subi}.html_fname=[opt.dirname '/plots/session' num2str(subnum) '/report.html'];    
        
    if ~exist(opt_report.sub_reports{subi}.html_fname,'file')
        opt_report.sub_reports{subi}.html_fname=[];
        opt_report.sub_reports{subi}.title=['Session ' num2str(subnum) ' (NOT FOUND)'];
    else,
        opt_report.sub_reports{subi}.title=['Session ' num2str(subnum) ];
    end;
end;

%%%%%%%%%%%%%%%%%%%%
%% diagnostic plots over all sessions
[opt opt_report]=opt_report_summary_plots(opt, opt_report);

%%%%%%%%%%%%%%%%%%%%
%% generate web report
opt.results.report=osl_report_write(opt_report);        

opt.fname=[opt.dirname '/opt'];

save(opt.fname, 'opt');
 
disp(['To view OPT report, point your browser to <a href="' opt.results.report.html_fname '">' opt.results.report.html_fname '</a>']);

diary off;
