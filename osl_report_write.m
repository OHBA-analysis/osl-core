function report=osl_report_write(report)
 
report.html_fname=[report.dir '/report.html'];
    
html=['<p>' report.title];

if(isfield(report,'logfile'))
    plotname = report.logfile(length(report.dir)+2:end);
    html=[html '<p><a href="' plotname '">Log file</a>'];
end;

% add links to sub reports:
if ~isempty(report.sub_reports)
    for rr=1:length(report.sub_reports),
        plotname = report.sub_reports{rr}.html_fname(length(report.dir)+2:end);
        html=[html '</p><p><a href="' [plotname] '">' report.sub_reports{rr}.title ' report</a>'];    
    end;
end;

% add any plots
if ~isempty(report.plot_names)
    for ii=1:length(report.plot_names),
        plotname = report.plot_names{ii}(length(report.dir)+2:end);
        html=[html '</p><p>' report.plot_titles{ii} '<p><img src="' plotname '.' report.fig_format '" width=650 >'];
    end;
end;

fid=fopen(report.html_fname,'w');
fprintf(fid,html);
fclose(fid);

end
