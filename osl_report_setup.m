function report=osl_report_setup(report_dir, title, logfile)
 
report=[];

if(nargin>2)
    report.logfile=logfile;
end

report.dir=report_dir;
report.title=title;

report.plot_fname=[];

report.index=0;    
report.plot_names={};
report.sub_reports={};

if ~exist(report_dir,'dir')
    mkdir(report_dir);
end

report.html_fname=[report.dir '/report.html'];

end
