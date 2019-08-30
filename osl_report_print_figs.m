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
    r = 150; % pixels per inch

    set(fighands(ii),'PaperPositionMode','auto');
    ppm = get(fighands(ii),'PaperPosition');
    set(fighands(ii),'PaperPosition',[0 0 ppm(3:4)]);
    set(fighands(ii),'PaperSize',ppm(3:4));

    if ~all(get(fighands(ii),'Color')==0.94) % The user set a custom background color
        set(fighands(ii),'InvertHardcopy','off'); % Use it
    else
        set(fighands(ii),'InvertHardcopy','on'); % Use it
    end

    report.plot_names{report.index}=[report.dir '/' fignames{ii}];
    warning off;
    print(fighands(ii), ['-d' report.fig_format], sprintf('-r%d',r), report.plot_names{report.index});
    warning on;
    %exportfig(fighands(ii),report.plots{report.index},'Color','rgb','Format',report.fig_format);
 
    report.plot_titles{report.index}=figtitles{ii};
    
    close(fighands(ii));        

end
