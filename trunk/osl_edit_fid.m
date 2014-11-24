function D = osl_edit_fid(D)
% Interactive tool to remove rogue headshape points stored in D.fiducials
% Adam Baker 2014

if ~isa(D,'meeg')
    try
        D = spm_eeg_load(D);
    catch
        error('Not a valid SPM file');
    end
end

fid = D.fiducials;

hf = figure('Toolbar','figure','Menubar','none','CloseRequestFcn',@closefig,'color','w','numbertitle','off','name','click or shift-click points to select, close figure to remove selected');
ha = axes('parent',hf);
plotfid
h = brush;
set(h,'Color',[1 0 0],'Enable','on');



    function closefig(src,evnt)
        selection = questdlg('Remove selected points?','Exit?');
        switch selection,
            case 'Yes'
                brushedData = get(findall(gca,'tag','Brushing'),'Xdata');
                brushedData = isnan(brushedData);
                fid.pnt(~brushedData,:) = [];
                D = fiducials(D,fid);
                plotfid
                pause(1)
                D.save;
            case 'No'
                
            case 'Cancel'
                return
        end
        delete(hf)
    end


    function plotfid
        plot3(fid.pnt(:,1),fid.pnt(:,2),fid.pnt(:,3),'color',[0.2 0.8 0.2],'marker','.','LineStyle','none','parent',ha)
        axis(ha,'image');
        axis(ha,'off');
    end

waitfor(hf)
D = spm_eeg_load(D);


end