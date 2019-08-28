function report=osl_report_add_text(report, txt, index)
    
    % add text to report at specified index in list of items
    % if no index provided then append on the end
    
    report.index=report.index+1;
    
    if nargin<3
        % append to the end
        index=report.index;
    else
        % insert
        titles={};
        if index>1
            for ii=1:index-1
                titles{ii}=report.plot_titles{ii};
            end
        end
        
        for ii=index+1:report.index
            titles{ii}=report.plot_titles{ii-1};
        end
        report.plot_titles=titles;
        
        names={};
        if index>1
            for ii=1:index-1
                names{ii}=report.plot_names{ii};
            end
        end
        
        for ii=index+1:report.index
            names{ii}=report.plot_names{ii-1};
        end
        report.plot_names=names;
    end
       
    report.plot_titles{index}=txt;
    report.plot_names{index}=[];

end
