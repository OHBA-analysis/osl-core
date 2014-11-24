function report=osl_report_setup(report_dir, report_title, report_logfile)
 
report=[];

if(nargin>2)
    report.logfile=report_logfile;
end;

report.dir=report_dir;
report.title=report_title;

report.plot_fname=[];

report.index=0;    
report.plot_names={};
report.sub_reports={};

warning off;
mkdir(report_dir);
warning on;

end
