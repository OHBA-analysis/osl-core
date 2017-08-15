function report=osl_report_print_figs(report)
 
if isstr(report.fig_name),
    fighands(1)=report.fig_handle;
    fignames{1}=report.fig_name;
    figtitles{1}=report.fig_title;
    
else
    fighands=report.fig_handle;
    fignames=report.fig_name;
    figtitles=report.fig_title;
end;

report.fig_format='png';

for ii=1:length(fighands),
    
    report.index=report.index+1;
    sfigure(fighands(ii));
    r = 150; % pixels per inch
    set(fighands(ii), 'PaperUnits', 'inches', 'PaperPosition', 2*get(fighands(ii),'Position')/r);
            
    report.plot_names{report.index}=[report.dir '/' fignames{ii}];
    print(fighands(ii), ['-d' report.fig_format], sprintf('-r%d',r), report.plot_names{report.index});
    %exportfig(fighands(ii),report.plots{report.index},'Color','rgb','Format',report.fig_format);
 
    report.plot_titles{report.index}=figtitles{ii};
    
    %close(fighands(ii));        

end
