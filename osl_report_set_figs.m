function report=osl_report_set_figs(report, fig_names, fig_handles, fig_titles)
 
    if nargin>2
        report.fig_handle=fig_handles;
        report.fig_name=fig_names;    
    else
        report.fig_handle=sfigure;
        report.fig_name=fig_names;    
    end
    
    if nargin>3
        report.fig_title=fig_titles;
    else
        report.fig_title=fig_names;
    end
end

