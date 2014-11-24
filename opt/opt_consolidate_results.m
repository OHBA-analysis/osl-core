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
opt.results.spm_files=cell(num_sessions,1);
opt.results.spm_files_basenames=cell(num_sessions,1);
opt.results.spm_files_epoched=cell(num_sessions,1);
opt.results.spm_files_epoched_basenames=cell(num_sessions,1);

for subi=1:length(opt.sessions_to_do), ii=opt.sessions_to_do(subi);
         
    % try and load in results
    opt.results.fnames{ii}=['session' num2str(ii)];
    
    try,
        opt_results=opt_load_results(opt, opt.results.fnames{ii});
    catch,
        opt.results.fnames{ii}=[];
    end;
    
    % build logfile
    runcmd(['cat ' opt_results.logfile '>' opt.results.logfile]);
    
    opt.results.spm_files_basenames{ii}=opt_results.spm_files_basename;
    opt.results.spm_files{ii}=[opt.dirname '/'  opt.results.spm_files_basenames{ii}];
 
    if opt.epoch.do,
    
        opt.results.spm_files_epoched_basenames{ii}=opt_results.spm_files_epoched_basename;
        opt.results.spm_files_epoched{ii}=[opt.dirname '/'  opt.results.spm_files_epoched_basenames{ii}];

    end;
end

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
