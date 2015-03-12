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

hf = figure('NumberTitle',                'off'     ,...
    'Menubar',               'none'    ,...
    'DockControls',          'off'     ,...
    'Color',                 'w'       ,...
    'CloseRequestFcn',       @closefig ,...
    'KeyPressFcn',           @key_press);

uitools.toolbar = uitoolbar;

[icon,~,a] = imread(fullfile(matlabroot,'toolbox','matlab','icons','file_save.png'));
a = double(a); a(a == 0) = nan; icon = double(icon) ./ 65535; icon(:,:,1) = icon(:,:,1) .* a ./ 65535;
uitools.save   = uipushtool(uitools.toolbar,'ClickedCallback',@save_meeg,'CData',icon,  'TooltipString','Save meeg object');

[icon,~,a]  = imread(fullfile(matlabroot,'toolbox','matlab','icons','tool_rotate_3d.png'));
a = double(a); a(a == 0) = nan; icon = double(icon) ./ 65535; icon(:,:,1) = icon(:,:,1) .* a ./ 65535;
uitools.rotate = uitoggletool(uitools.toolbar,'ClickedCallback',@togglerot,'CData',icon,'TooltipString','Rotation');

[icon,~,a]   = imread(fullfile(matlabroot,'toolbox','matlab','icons','tool_data_brush.png'));
a = double(a); a(a == 0) = nan; icon = double(icon) ./ 65535; icon(:,:,1) = icon(:,:,1) .* a ./ 65535;
uitools.brush  = uitoggletool(uitools.toolbar,'ClickedCallback',@togglebrush,'CData',icon,'TooltipString','Select data');


ha = axes('parent',hf);
plotfid
h = brush;
set(h,'Color',[1 0 0],'Enable','on');

    function closefig(~,~)
        delete(hf)
    end

    function save_meeg(~,~)
        if any(get(findall(gca,'tag','Brushing'),'Xdata'))
            selection = questdlg('Remove selected points?','Yes','No');
            switch selection,
                case 'Yes'
                    brushedData = get(findall(gca,'tag','Brushing'),'Xdata');
                    brushedData = isnan(brushedData);
                    fid.pnt(~brushedData,:) = [];
                    D = fiducials(D,fid);
                    [AZ,EL] = view;
                    plotfid
                    view(AZ,EL);
                    pause(1)
                    D.save;
                case 'No'
                    
                case 'Cancel'
                    return
            end
        end
    end


    function togglerot(~,~)
        brush off
        set(uitools.brush,'state','off')
        rotate3d
    end

    function togglebrush(~,~)
        rotate3d off
        set(uitools.rotate,'state','off')
        brush
    end

    function plotfid
        plot3(fid.pnt(:,1),fid.pnt(:,2),fid.pnt(:,3),'color',[0.2 0.8 0.2],'marker','.','LineStyle','none','parent',ha)
        axis(ha,'image');
        axis(ha,'off');
    end


end
