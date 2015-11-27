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
hfid = [];

hf = figure('NumberTitle',                'off'     ,...
    'Menubar',               'none'    ,...
    'DockControls',          'off'     ,...
    'Color',                 'k'       ,...
    'CloseRequestFcn',       @closefig ,...
    'KeyPressFcn',           @key_press,...
    'Name',                  fullfile(D.path,D.fname));

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
brushedData = [];
hBrush = brush;

togglerot;
set(uitools.rotate,'state','on');

    function closefig(~,~)
        delete(hf)
    end

    function save_meeg(~,~)
        brushedData = get_brushed_data;
        if any(brushedData)
            selection = questdlg('Remove selected points?');
            switch selection,
                case 'Yes'
                    fid.pnt(brushedData,:) = [];                  
                    if isfield(fid,'label')
                        fid.label(brushedData) = [];
                    end
                    D = fiducials(D,fid);
                    [AZ,EL] = view;
                    plotfid
                    view(AZ,EL);
                    pause(0.1)
                    D.save
                    D = spm_eeg_load(fullfile(D.path,D.fname));
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
        hfid = plot3(fid.pnt(:,1),fid.pnt(:,2),fid.pnt(:,3),'color','g','marker','.','LineStyle','none','parent',ha);
        axis(ha,'image');
        axis(ha,'off');
    end

    function brushedData = get_brushed_data
        
        if verLessThan('matlab', '8.4') % is HG1
            brushedData = get(findall(gca,'tag','Brushing'),'Xdata');
            brushedData = ~isnan(brushedData);
        else % is HG2
            brushedData = hfid.BrushData;
        end
            
    end
end
















