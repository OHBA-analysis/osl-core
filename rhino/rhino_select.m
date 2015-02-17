function fid_positions = rhino_select(mesh,fid_positions,fid_labels)

MainFig = figure('NumberTitle',           'off'     ,...
                 'Menubar',               'none'    ,...
                 'DockControls',          'off'     ,...
                 'CloseRequestFcn',       @closefig ,...
                 'KeyPressFcn',           @key_press);
                       
uitools.toolbar = uitoolbar;
rotateicon = load(fullfile(matlabroot,'toolbox','matlab','icons','rotate.mat'));
uitools.rotate = uitoggletool(uitools.toolbar,'ClickedCallback',@togglerot,'CData',rotateicon.cdata,'TooltipString','Toggle rotation');

if any(strcmp(get(get(MainFig,'children'),'type'),'axes'))
    [az,el] = view;
else
    az = 90; el = 0;
end

ha = axes('parent',MainFig); hold on


set(MainFig,'WindowButtonMotionFcn',@movelight,'Interruptible','off','busyaction','cancel');

rotate3d(MainFig);

hPatch = patch(struct('faces',mesh.faces,'vertices',mesh.vertices),'FaceColor',[238,206,179]./255,'EdgeColor','none','Parent',ha);

view([az,el]);

axis(ha,'image','off')
material dull
lighting gouraud
hl = camlight('headlight');

% Note to anyone who sees this code - it was a massive ball ache to get
% a single update function for each datatip, and in the end I couldn't
% be arsed so I've just made a separate one for each fiducial. If you
% can find a nicer solution then you're a better man than I am!
cursorNAS = datacursormode(MainFig);
set(cursorNAS, 'enable','on', 'SnapToDataVertex','on', 'UpdateFcn',@getFidPositionNAS, 'NewDataCursorOnClick',false);
hcursorNAS = cursorNAS.createDatatip(handle(hPatch));
set(hcursorNAS,'HandleVisibility','off','UIContextMenu',[]);

cursorLPA = datacursormode(MainFig);
set(cursorLPA, 'enable','on', 'SnapToDataVertex','on', 'UpdateFcn',@getFidPositionLPA, 'NewDataCursorOnClick',false);
hcursorLPA = cursorLPA.createDatatip(handle(hPatch));
set(hcursorLPA,'HandleVisibility','off','UIContextMenu',[]);

cursorRPA = datacursormode(MainFig);
set(cursorRPA, 'enable','on', 'SnapToDataVertex','on', 'UpdateFcn',@getFidPositionRPA, 'NewDataCursorOnClick',false);
hcursorRPA = cursorRPA.createDatatip(handle(hPatch));
set(hcursorRPA,'HandleVisibility','off','UIContextMenu',[]);

[~,pointNAS] = min(sum((mesh.vertices - repmat(fid_positions(strncmpi(fid_labels,'nas',3),:),size(mesh.vertices,1),1)).^2,2));
[~,pointLPA] = min(sum((mesh.vertices - repmat(fid_positions(strncmpi(fid_labels,'lpa',3),:),size(mesh.vertices,1),1)).^2,2));
[~,pointRPA] = min(sum((mesh.vertices - repmat(fid_positions(strncmpi(fid_labels,'rpa',3),:),size(mesh.vertices,1),1)).^2,2));

%min(sum((mesh.vertices - repmat(fid_positions(strncmpi(fid_labels,'nas',3),:),size(mesh.vertices,1),1)).^2,2))

set(get(hcursorNAS,'DataCursor'),'Position',fid_positions(strncmpi(fid_labels,'nas',3),:));
set(get(hcursorLPA,'DataCursor'),'Position',fid_positions(strncmpi(fid_labels,'lpa',3),:));
set(get(hcursorRPA,'DataCursor'),'Position',fid_positions(strncmpi(fid_labels,'rpa',3),:));

set(get(hcursorNAS,'DataCursor'),'DataIndex',pointNAS,'TargetPoint',mesh.vertices(pointNAS,:));
set(get(hcursorLPA,'DataCursor'),'DataIndex',pointLPA,'TargetPoint',mesh.vertices(pointLPA,:));
set(get(hcursorRPA,'DataCursor'),'DataIndex',pointRPA,'TargetPoint',mesh.vertices(pointRPA,:));

uiwait(MainFig)


    function movelight(varargin)
        delete(hl)
        hl = camlight('headlight');
    end

    function txt = getFidPositionNAS(~,evnt)           
        pos = get(evnt,'Position');
        lbl = char(fid_labels(strncmpi(fid_labels,'nas',3)));
        txt =  {lbl,     ...
               ['x: ',num2str(pos(1))], ...
               ['y: ',num2str(pos(2))], ...
               ['z: ',num2str(pos(3))]}; 
    end

    function txt = getFidPositionLPA(~,evnt)           
        pos = get(evnt,'Position');
        lbl = char(fid_labels(strncmpi(fid_labels,'lpa',3)));
        txt =  {lbl,     ...
               ['x: ',num2str(pos(1))], ...
               ['y: ',num2str(pos(2))], ...
               ['z: ',num2str(pos(3))]}; 
    end

    function txt = getFidPositionRPA(~,evnt)           
        pos = get(evnt,'Position');
        lbl = char(fid_labels(strncmpi(fid_labels,'rpa',3)));
        txt =  {lbl,     ...
               ['x: ',num2str(pos(1))], ...
               ['y: ',num2str(pos(2))], ...
               ['z: ',num2str(pos(3))]}; 
    end


    function togglerot(varargin)
        rotate3d
    end

    function closefig(varargin)

        try
            cursors = datacursormode(MainFig);
            cursorLabels = arrayfun(@(c) cursors.DataCursors(c).String{1},1:3,'UniformOutput',0);
            
            for lbl = fid_labels
                fid_positions(strcmp(lbl,fid_labels),:) = get(get(cursors.DataCursors(strcmp(cursorLabels,lbl)), ...
                    'DataCursor'),'Position');
            end
            
            delete(MainFig)
            
        catch
            delete(MainFig)
            error('error while closing figure')
        end
    end


end
