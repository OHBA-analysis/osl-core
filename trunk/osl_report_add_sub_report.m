function report=osl_report_add_sub_report(report, sub_report)
 
num=length(report.sub_reports);
report.sub_reports{num+1}=sub_report;
    
end