function D = oslview(D)
    % MEG continuous data viewer
    % OSLVIEW(D)
    % ----------------------------------------------------------------
    % D - SPM meeg object or file name
    % ----------------------------------------------------------------
    %
    % RA 2017
    % Updated 21/10/14
    % AB 2012
    %
    %
    % TODO:
    % - go to channel
    % - epoched data?
    % - Have a reset view button or right click option always available.
    % - Tighter control over SideWindow bars

    if nargout == 0
        fprintf(2,'Warning - oslview returns an MEEG object but it is not being assigned to a variable\n');
    end

    assert(strcmp(D.type,'continuous'),'oslview is only compatible with continuous MEEG objects');

    % Get directory of the viewer
    viewer_dir = strrep(which('oslview'),'oslview.m','');

    % Initialise shared variables
    drag = struct('state',0); % Struct to track the state of the main window drag
    chan_time = []; % Array of time points for plots on main window
    chan_data = []; % Matrix of value points for plots on main window
    chan_inds = []; % Array specifying order in which to display the components
    chan_labels = D.chanlabels; % Cache this because accessing it dynamically is extremely slow

    selection = struct('line',[],'bar',[]); % Temporary lines for context menu

    Nchannels = [];
    SideWindowData  = [];
    PanWindowData = [];
    offsets = [];
    yzoom = [0 1];
    BadEpochs = {};
    chanbar = []; % Handles to bars in sidebar
    chansig = []; % Handle to line plots in main window
    badevents_line = [];
    badevents_patch = [];
    PanWindow_line = [];
    PanWindow_box = [];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create UI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MainFig = figure('Name',['OSLview - ' strtok(D.fname,'.')],...
                        'NumberTitle','off',...
                        'Menubar','none',...
                        'DockControls','off',...
                        'WindowButtonDownFcn',@btn_down ,...
                        'WindowButtonUpFcn',@btn_up ,...
                        'KeyPressFcn',@key_press,...
                        'ResizeFcn',@(a,b,c) createlayout,...
                        'Units','pixels',...
                        'Visible','off');
            
    c = onCleanup(@() delete(MainFig));
                
    % Create plotting windows
    MainWindow = axes('parent',MainFig, 'units','pixels','XGrid','on');
    PanWindow   = axes('parent',MainFig, 'units','pixels');
    SideWindow = axes('parent',MainFig, 'units','pixels');

    % Load tool icons
    icons = load([viewer_dir 'oslview.mat']);

    % Create toolbars
    uitools.toolbar     = uitoolbar;
    uitools.expand          = uipushtool(uitools.toolbar,       'ClickedCallback',@inc_xwidth,          'CData',icons.expand,                   'TooltipString','Increase time window');
    uitools.shrink          = uipushtool(uitools.toolbar,       'ClickedCallback',@dec_xwidth,          'CData',icons.shrink,                   'TooltipString','Decrease time window');
    uitools.rect                = uitoggletool(uitools.toolbar, 'ClickedCallback',@mouse_mode,          'CData',icons.rect,                     'TooltipString','Zoom to channels');
    uitools.zoomin          = uipushtool(uitools.toolbar,       'ClickedCallback',@inc_scale,               'CData',icons.zoomin,                   'TooltipString','Increase scale');
    uitools.zoomout     = uipushtool(uitools.toolbar,       'ClickedCallback',@dec_scale,               'CData',icons.zoomout,              'TooltipString','Decrease scale');

    % Create menu for channel selection
    uitools.menu_channels = uimenu('label','Channels');
    channel_types = unique(D.chantype);
    for chtype = 1:length(channel_types)
        uimenu(uitools.menu_channels,'label',channel_types{chtype},'Callback',@switch_chantype);
        [BadEpochs.(channel_types{chtype}),InheritedBadEpochs.(channel_types{chtype})] = read_bad_events(D,channel_types{chtype});
    end
        

    % Create Context Menu - Main Window
    ContextMenuM.Menu                   = uicontextmenu('callback',@ContextMenuOn);
    ContextMenuM.ChannelLabel = uimenu(ContextMenuM.Menu, 'label','No Channel Selected','Enable','off');
    ContextMenuM.MarkEvent      = uimenu(ContextMenuM.Menu, 'label','Mark Event'                , 'callback',@mark_bad);
    ContextMenuM.SetBad             = uimenu(ContextMenuM.Menu, 'label','Set Channel as Bad', 'callback',@bad_channel);
    set(MainWindow,'uicontextmenu',ContextMenuM.Menu);

    % Create Context Menu - Side Window axis
    ContextMenuS.Menu           = uicontextmenu;
    ContextMenuS.Reorder    = uimenu(ContextMenuS.Menu, 'label','Sort by metric','callback',@cb_ContextMenuS);
    ContextMenuS.Variance = uimenu(ContextMenuS.Menu, 'label','Variance','Checked','on','callback',@cb_ContextMenuS);
    ContextMenuS.Kurtosis = uimenu(ContextMenuS.Menu, 'label','Kurtosis','Checked','off','callback',@cb_ContextMenuS);
    set(SideWindow,'uicontextmenu',ContextMenuS.Menu);

    % Create Context Menu - Side Window patch
    ContextMenuSP.Menu = uicontextmenu('callback',@ContextMenuOnSP);
    ContextMenuSP.ChannelLabel = uimenu(ContextMenuSP.Menu, 'label','No Channel Selected','Enable','off');
    ContextMenuSP.SetBad = uimenu(ContextMenuSP.Menu, 'label','Set Channel as Bad', 'callback',@bad_channel);

    % Setup figure layout and plot data
    createlayout

    % Get sample times
    Nsamples = D.nsamples;
    t = D.time;

    % Set gain for zooming
    G = 1;

    % Set zooming factor
    Ginc = 2;

    % Handle for zoom control
    hz = zoom(MainFig);

    % Set initial channel type
    if any(strcmp(D.chantype,'MEGMAG')) && any(strcmp(D.chantype,'MEGPLANAR')) % ELEKTA
        channel_type = 'MEGPLANAR';
    elseif any(strcmp(D.chantype,'MEGMAG')) % 4D
        channel_type = 'MEGMAG';
    elseif any(strcmp(D.chantype,'MEGGRAD')) % CTF
        channel_type = 'MEGGRAD';
    elseif any(strcmp(D.chantype,'MEG'))
        channel_type = 'MEG';
    else % catch unrecognised case
        channel_type = D.chantype{1};
    end
    channel_setup; % will also call redraw & redraw_Sidewindow

    set(MainFig,'visible','on');

    drawnow
    warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
    jFig = get(handle(MainFig),'JavaFrame');
    jFig.setMaximized(true);
    warning on MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame

    %drag_listener1 = addlistener(MainWindow,'XLim','PostSet',@redraw);
    set(PanWindow,'ButtonDownFcn',@pan_jump)

    uiwait(MainFig) % Wait for user to close window

    for chtype = 1:length(channel_types)
        D = set_bad_events(D,BadEpochs.(channel_types{chtype}),channel_types{chtype});
    end

    return

    function createlayout(varargin)
        % Perform one-off creation of axes objects etc.
        Layout.MainFig.position = get(MainFig,'position');
        Layout.FigWidth = Layout.MainFig.position(3);
        Layout.FigHeight = Layout.MainFig.position(4);
        Layout.PanWindowWidth = 50;
        Layout.SideWindowWidth = 50;
        Layout.BorderWidth = 10;
        Layout.LabelWidth = 20;
        
        
        % Main Window
        Layout.MainWindow.X(1) = Layout.BorderWidth + Layout.LabelWidth;
        Layout.MainWindow.X(2) = Layout.FigWidth - 2*Layout.BorderWidth - Layout.SideWindowWidth;
        Layout.MainWindow.Y(1) = 2*Layout.BorderWidth + Layout.PanWindowWidth;
        Layout.MainWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth;
        Layout.MainWindow.position = [Layout.MainWindow.X(1),Layout.MainWindow.Y(1) ,abs(diff(Layout.MainWindow.X)),abs(diff(Layout.MainWindow.Y))];

        
        % Pan Window
        Layout.PanWindow.X(1) = Layout.BorderWidth + Layout.LabelWidth;
        Layout.PanWindow.X(2) = Layout.FigWidth - 2*Layout.BorderWidth - Layout.SideWindowWidth;
        Layout.PanWindow.Y(1) = Layout.BorderWidth;
        Layout.PanWindow.Y(2) = Layout.PanWindowWidth;
        Layout.PanWindow.position = [Layout.PanWindow.X(1),Layout.PanWindow.Y(1),abs(diff(Layout.PanWindow.X)),abs(diff(Layout.PanWindow.Y))];
        
            
        % Side Window
        Layout.SideWindow.X(1) = Layout.FigWidth - Layout.BorderWidth - Layout.SideWindowWidth;
        Layout.SideWindow.X(2) = Layout.FigWidth - Layout.BorderWidth;
        Layout.SideWindow.Y(1) = 2*Layout.BorderWidth + Layout.PanWindowWidth;
        Layout.SideWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth;
        Layout.SideWindow.position = [Layout.SideWindow.X(1),Layout.SideWindow.Y(1),abs(diff(Layout.SideWindow.X)),abs(diff(Layout.SideWindow.Y))];

        set(MainWindow,'position',Layout.MainWindow.position,'XLim',D.time([1 end]));
        set(PanWindow,'position',Layout.PanWindow.position);
        set(SideWindow,'position',Layout.SideWindow.position);
        box(MainWindow,'on')
        box(PanWindow,'on')
        box(SideWindow,'on')

    end

    function pan_jump(varargin)
        % Reposition the window by clicking on bottom plot
        xspan = get(MainWindow,'XLim');
        set(MainWindow,'XLim',varargin{2}.IntersectionPoint(1)+[-0.5 0.5]*diff(xspan));
        redraw
    end

    function xspan = validate_xspan(xspan)
        % If we want to set the axis limits, it needs to be clipped to the data
        xwidth = xspan(2)-xspan(1);
        if xwidth >= (D.time(end)-D.time(1))
            xspan = [D.time(1) D.time(end)];
            return
        end
        
        if xspan(1) < D.time(1)
            xspan = xspan - xspan(1) + D.time(1);
        end

        if xspan(end) > D.time(end)
            xspan = xspan - (xspan(end)-D.time(end));
        end

    end

    function redraw
        % Update all plots

        ContextMenuOff() % If redrawing, then we need to get rid of the context menu temporary lines
        
        xspan = validate_xspan(get(MainWindow,'XLim'));
        set(MainWindow,'XLim',xspan);

        % Based on the time window desired, calculate the indexes of the D object to sample
        x_window = round(xspan*D.fsample - D.time(1)); % Convert from time to samples
        x_window(1) = max(x_window(1),1);
        x_window(2) = min(x_window(2),length(D.time));
        xs = unique(round(linspace(x_window(1),x_window(2),1000)));
        
        ylim_mainwindow = diff(get_ylims).*yzoom + min(get_ylims);
        
        % Plot channel signals within current range
        chan_time = t(xs);
        chan_data = G*ones(size(Nchannels)).*D(chan_inds,xs) + repmat(offsets(:),size(xs));

        for ch = 1:length(chansig)
            set(chansig(ch),'XData',chan_time,'YData',chan_data(ch,:),'LineWidth',0.5,'LineStyle','-','tag',chan_labels{chan_inds(ch)});
        end

        set(MainWindow,'ylim',ylim_mainwindow)
        set(MainWindow,'ytick',[],'yticklabel',[])

        % Add context menu to plotted signals and label with channel name
        set(chansig,'uicontextmenu',ContextMenuM.Menu);
        
        % Replace bad channels with dashed lines
        ch_bad = get_bad_channels;
        
        set(chansig(ch_bad),'linestyle','--','linewidth',1);

        %% PAN WINDOW
        set(PanWindow_line,'XData',t,'YData',PanWindowData);
        if max(PanWindowData) == 0
            ylim = [0 1];
        else
            ylim = [0 max(PanWindowData)];
        end

        set(PanWindow,'YLim',ylim);
        
        box_x = [t(xs(1)) t(xs(end)) t(xs(end)) t(xs(1))];
        box_y = [ylim(1) ylim(1) ylim(2) ylim(2)];
        set(PanWindow_box,'XData',box_x,'YData',box_y);

        %% SIDE WINDOW
        
        % Read offsets for the bar positions from elsewhere as well
        ch_bad = get_bad_channels;
        SideWindowData_plot = SideWindowData;
        SideWindowData_plot(ch_bad) = NaN;

        % Because offsets are set to cumsum when they are computed, they are guaranteed to be sorted
        % So they can be replaced here with just simple offsets
        sidebar_offsets = (1:length(offsets))-0.5;
        for ch = 1:Nchannels
            bar_patch(chanbar(ch),sidebar_offsets(ch),SideWindowData_plot(ch),1,chan_labels{chan_inds(ch)});
        end
        
        set(SideWindow,'ylim',[0 length(offsets)])
        if length(SideWindowData_plot) == 1
            set(SideWindow,'xlim',[0 1])
        else
            set(SideWindow,'xlim',[min(SideWindowData_plot) max(SideWindowData_plot)])
        end
        set(SideWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[]);

        if strcmp(get(ContextMenuS.Variance,'Checked'),'on')
            xlabel(SideWindow,'Variance');
        else
            xlabel(SideWindow,'Kurtosis');
        end

        %% BAD EVENTS
        yls= [get(MainWindow,'YLim') get(PanWindow,'YLim')];
        yl = [min(yls) max(yls)];

        green_x = NaN;
        green_y = NaN;
        red_x = NaN;
        red_y = NaN;
        patch_x = [];
        patch_y = [];

        for j = 1:length(BadEpochs.(channel_type))
            green_x = [green_x BadEpochs.(channel_type){j}(1) BadEpochs.(channel_type){j}(1) NaN];
            green_y = [green_y yl NaN];
            if numel(BadEpochs.(channel_type){j}) == 2
                red_x = [red_x BadEpochs.(channel_type){j}(2) BadEpochs.(channel_type){j}(2) NaN];
                red_y = [red_y yl NaN];
                patch_x(:,end+1) = [BadEpochs.(channel_type){j}(1);BadEpochs.(channel_type){j}(2);BadEpochs.(channel_type){j}(2);BadEpochs.(channel_type){j}(1)];
                patch_y(:,end+1) = [yl(1);yl(1);yl(2);yl(2)];
            end
        end

        for j = 1:length(InheritedBadEpochs.(channel_type))
            if numel(InheritedBadEpochs.(channel_type){j}) == 2
                patch_x(:,end+1) = [InheritedBadEpochs.(channel_type){j}(1);InheritedBadEpochs.(channel_type){j}(2);InheritedBadEpochs.(channel_type){j}(2);InheritedBadEpochs.(channel_type){j}(1)];
                patch_y(:,end+1) = [yl(1);yl(1);yl(2);yl(2)];
            end
        end

        set(badevents_line([1,3]),'XData',green_x,'YData',green_y);
        set(badevents_line([2,4]),'XData',red_x,'YData',red_y);
        set(badevents_patch,'XData',patch_x,'YData',patch_y);
        
    end

    function bar_patch(h_patch,offset,value,width,tag)
        % Helper function to use a patch object as a fast horizontal bar
        xv = [0 value];
        yv = offset + [-1 1]*width/2;
        set(h_patch,'YData',[yv(1) yv(1) yv(2) yv(2)],'XData',[xv(1) xv(2) xv(2) xv(1)],'tag',tag)
    end

    function switch_chantype(src,~)
        channel_type = get(src,'label');

        % Update bad segments so that inherited artefacts are correctly displayed
        for chtype = 1:length(channel_types)
            D = set_bad_events(D,BadEpochs.(channel_types{chtype}),channel_types{chtype});
        end     
        for chtype = 1:length(channel_types)
            [BadEpochs.(channel_types{chtype}),InheritedBadEpochs.(channel_types{chtype})] = read_bad_events(D,channel_types{chtype});
        end

        % Update the display
        channel_setup;
    end

    function channel_setup
        % This function clears the axes and makes new line, bar, and patch objects
        % Occurs when the channel type is changed
        % redraw() edits the properties of the handles created here

        % Update the channel statistics
        chan_inds = find(strcmp(D.chantype, channel_type));

        calcPlotStats

        if strcmp(get(ContextMenuS.Reorder,'Checked'),'on')
            [~,idx] = sort(SideWindowData,'ascend');
            chan_inds = chan_inds(idx);
        else
            chan_inds = sort(chan_inds,'ascend');
        end

        Nchannels = length(chan_inds);

        % Set channel labels & plot colours
        chancols = colormap(lines); 
        chancols = chancols(1:7,:);
        chancols = repmat(chancols,ceil(Nchannels/7),1);
        
        % Calculate side and pan window statistics
        SideWindowData = [];
        PanWindowData = [];
        calcPlotStats
        
        % Make all of the bars
        xspan = get(MainWindow,'XLim');
        cla(MainWindow);
        hold(MainWindow,'on')
        chansig = plot(MainWindow,xspan,ones(2,Nchannels));
        badevents_line(1) = plot(MainWindow,[NaN NaN],[NaN NaN],'g','LineWidth',2,'LineStyle','--','HitTest','off','YLimInclude','off');
        badevents_line(2) = plot(MainWindow,[NaN NaN],[NaN NaN],'r','LineWidth',2,'LineStyle','--','HitTest','off','YLimInclude','off');
        badevents_patch(1) = patch(nan(4,1),nan(4,1),'k','LineStyle','none','Parent',MainWindow,'FaceAlpha',0.1,'HitTest','off','YLimInclude','off');

        cla(SideWindow);
        hold(SideWindow,'on')
        for ch = 1:Nchannels
            chanbar(ch) = patch(zeros(1,4),zeros(1,4),chancols(ch,:),'Parent',SideWindow,'uicontextmenu',ContextMenuSP.Menu,'LineStyle','none');
        end

        cla(PanWindow);
        hold(PanWindow,'on');
        set(PanWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[],'ButtonDownFcn',@pan_jump,'XLim',[D.time(1) D.time(end)]);
        PanWindow_line = plot(PanWindow,NaN,NaN,'b','HitTest','off');
        PanWindow_box = patch(NaN,NaN,'r','parent',PanWindow,'facealpha',0.3,'HitTest','off','YLimInclude','off');
        badevents_line(3) = plot(PanWindow,[NaN NaN],[NaN NaN],'g','LineWidth',1,'LineStyle','-','HitTest','off','YLimInclude','off');
        badevents_line(4) = plot(PanWindow,[NaN NaN],[NaN NaN],'r','LineWidth',1,'LineStyle','-','HitTest','off','YLimInclude','off');
        badevents_patch(2) = patch(nan(4,1),nan(4,1),'k','Parent',PanWindow,'LineStyle','none','FaceAlpha',0.1,'HitTest','off','YLimInclude','off');

        redraw

    end

    function ContextMenuOn(varargin)
        
        cp = get(MainWindow,'CurrentPoint');
        
        % Find the time index of clicked time
        t_idx = find(chan_time<=cp(1),1,'last');

        % Find the channel with the closest value at that time 
        [~,ch] = min(abs(cp(1,2)-chan_data(:,t_idx)));
    
        set(ContextMenuM.ChannelLabel,'label',chan_labels{chan_inds(ch)});
        
        % Change context menu depending on whether channel is marked as bad
        if isempty(find(D.badchannels==find(strcmp(D.chanlabels,get(ContextMenuM.ChannelLabel,'label'))),1))
            set(ContextMenuM.SetBad,'label','Set Channel as Bad');
        else
            set(ContextMenuM.SetBad,'label','Set Channel as Good');
        end
        
        % Change context menu depending on whether time point with bad epoch
        if any(cellfun(@(x) ~sum(sign(x-cp(1))),BadEpochs.(channel_type)))
            set(ContextMenuM.MarkEvent,'label','Remove Event');
        else
            set(ContextMenuM.MarkEvent,'label','Mark Event');
        end

        if ishandle(selection.line)
            delete(selection.line);
        end
        selection.line = plot(MainWindow,chan_time,chan_data(ch,:),'Color','k','Linewidth',2);
    end


    function ContextMenuOnSP(varargin)
        ch = find(chanbar == gco); % This context menu is activated by clicking on a bar, so the patch should now be the current object
    
        set(ContextMenuSP.ChannelLabel,'label',chan_labels{chan_inds(ch)});
        
        % Change context menu depending on whether channel is marked as bad
        if isempty(find(D.badchannels==find(strcmp(D.chanlabels,get(ContextMenuSP.ChannelLabel,'label'))),1))
            set(ContextMenuSP.SetBad,'label','Set Channel as Bad');
        else
            set(ContextMenuSP.SetBad,'label','Set Channel as Good');
        end

        if ishandle(selection.line)
            delete(selection.line);
        end
        selection.line = plot(MainWindow,chan_time,chan_data(ch,:),'Color','k','Linewidth',2);
    end

    function ContextMenuOff(varargin)
        if ishandle(selection.line)
            delete(selection.line);
        end
    end

    function cb_ContextMenuS(src,~)

        switch src
            case ContextMenuS.Reorder
                if strcmp(get(ContextMenuS.Reorder,'Checked'),'off')
                    set(ContextMenuS.Reorder,'Checked','on');
                else
                    set(ContextMenuS.Reorder,'Checked','off');
                end

            case ContextMenuS.Variance
                set(ContextMenuS.Variance,'Checked','on');
                set(ContextMenuS.Kurtosis,'Checked','off');
            case ContextMenuS.Kurtosis
                set(ContextMenuS.Kurtosis,'Checked','on');
                set(ContextMenuS.Variance,'Checked','off');
        end

        channel_setup % Re-run channel setup because chan_inds might have changed
        
    end

    function btn_down(varargin)
        % Arm the drag mode by clicking on main plot

        % Did we click on the main figure?
        loc = get(MainFig,'CurrentPoint');
        pos = get(MainWindow,'Position');
        inside = loc(1) >= pos(1) && loc(1) <= pos(1)+pos(3) && loc(2) >= pos(2) && loc(2) <= pos(2)+pos(4);

        if ~inside
            drag.state = 0;
            return;
        end

        if strcmp(get(MainFig,'SelectionType'),'extend')
            mark_bad
            return
        end

        drag = struct;
        drag.initial_point = get(MainWindow,'CurrentPoint');
        drag.initial_xlim = get(MainWindow,'XLim');
        drag.state = 2; % 0 - no drag, 1 - drag in progress, 2 - drag armed

        ContextMenuOff % deselect highlighted channel
        set(MainFig,'WindowButtonMotionFcn',@btn_move)

        % yzoom = (get(MainWindow,'ylim') - min(get_ylims))./(diff(get_ylims));
    end

    function btn_move(varargin)
        if is_multiple_call
            return
        end

        p = get(MainWindow,'CurrentPoint');
        current_center = mean(get(MainWindow,'XLim')); % Current x limit
        initial_center = mean(drag.initial_xlim); % 
        displacement = p(1)-drag.initial_point(1)-(current_center-initial_center);

        % If drag is armed and needs to be initiated, set drag mode and switch to low res plot
        if drag.state == 2 && abs(displacement/diff(drag.initial_xlim))>0.025 % Start drag after moving just a little
            drag.state = 1;
            xs = unique(round(D.fsample*linspace(D.time(1),D.time(end),2000)));
            xs = xs(2:end);
            chan_data = G*ones(size(Nchannels)).*D(chan_inds,xs) + repmat(offsets(:),size(xs));
            for ch = 1:length(chansig)
                set(chansig(ch),'XData',t(xs),'YData',chan_data(ch,:),'LineWidth',0.5,'LineStyle','-','tag',chan_labels{chan_inds(ch)});
            end
        end

        if drag.state == 1
            new_center = initial_center-displacement;
            xspan = new_center + 0.5*[-1 1]*diff(drag.initial_xlim); 
            xspan = validate_xspan(xspan);
            set(MainWindow,'XLim',xspan);
            set(PanWindow_box,'XData',[xspan(1) xspan(2) xspan(2) xspan(1)]);
        end
    end

    function btn_up(varargin)
        % Disable drag movement when button is lifted
        set(MainFig,'WindowButtonMotionFcn',[]);
        if drag.state == 1
            redraw
        end
        drag.state = 0;
    end

    function mark_bad(varargin)
        % This callback adds event start/stop times and handles manipulation
        % of BadEpochs

        t_current = get(MainWindow,'CurrentPoint');
        t_current = t_current(1,1);
        t_window    = get(MainWindow,'xlim');

        if isempty(varargin)
            if any(cellfun(@(x) ~sum(sign(x-t_current(1))),BadEpochs.(channel_type))) && ~isempty(BadEpochs.(channel_type)) && length(BadEpochs.(channel_type){end}) ~= 1
                new_event = false;
            else
                new_event = true;
            end
        else
            if strcmp(get(ContextMenuM.MarkEvent,'label'),'Mark Event')
                new_event = true;
            else
                new_event = false;
            end
        end

        if new_event %  Create new event
            % if position of marker is within the first or last 1% of the window
            % then assume it was meant to be at the start or end.
            if t_current < t_window(1) + 0.01*diff(t_window)
                t_current = t_window(1);
            elseif t_current > t_window(1) + 0.99*diff(t_window)
                t_current = t_window(2);
            end

            if isempty(BadEpochs.(channel_type)) || length(BadEpochs.(channel_type){end}) == 2 % Start a new epoch
                BadEpochs.(channel_type){end+1}(1) = t_current;
                redraw
            else % End previous epoch, recompute stats
                BadEpochs.(channel_type){end}(2) = t_current;
                BadEpochs.(channel_type){end} = sort(BadEpochs.(channel_type){end});

                % Recompute the events
                t_bad = get_bad_inds();

                db = find(diff([0 t_bad]));
                onset = db(1:2:end);
                offset = db(2:2:end);
                if length(offset)<length(onset) 
                    offset(end+1) = length(t_bad); 
                end

                BadEpochs.(channel_type) = {};
                for j = 1:length(onset)
                    BadEpochs.(channel_type){j} = D.time([onset(j),offset(j)]);
                end
                
                calcPlotStats
                redraw
            end
        else    %   Remove existing event
            BadEpochs.(channel_type)(cellfun(@(x) ~sum(sign(x-t_current)),BadEpochs.(channel_type))) = [];
            calcPlotStats
            redraw
        end
    end


    function bad_channel(varargin)
        % Toggle the bad_channel state of the channel with the specified name

        % Read the bad channel from the parent 
        label = get(findobj(get(varargin{1},'Parent'),'Enable','off'),'Label'); % Get the parent ContextMenu, find the entry that's greyed out, and return its label

        bad_chan = find(strcmp(D.chanlabels,label));

        if isempty(find(D.badchannels==bad_chan,1))
            D = badchannels(D,bad_chan,1);
        else
            D = badchannels(D,bad_chan,0);
        end
            
        calcPlotStats()
        if ~isempty(selection.line)
            delete(selection.line);
        end
        redraw
    end


    function calcPlotStats()
        
        % Get statistic from side window context menu
        % Statistic is either variance or Kurtosis
        if strcmp(get(ContextMenuS.Variance,'Checked'),'on')
            stat = @nanvar;
        else
            stat = @kurtosis;
        end     
        
        % Bad Channel indices
        ch_bad = get_bad_channels;
        
        % Bad Sample indices
        t_bad = get_bad_inds(false) | get_bad_inds(true); % Exclude both direct and inherited times from the stats

        Dsig_inds = 1:20:D.nsamples; Dsig_inds(t_bad(Dsig_inds)) = [];
%       Dsig = std(D(:,Dsig_inds,:),[],2);
        offsets = 3*iqr(D(:,Dsig_inds,:),2);
        offsets(offsets==0) = eps;
        offsets(D.badchannels) = nan;
        offsets = offsets(chan_inds);
        offsets(isnan(offsets)) = nanmean(offsets);
        offsets = cumsum(offsets);
    
        % For the PanWindow, display the variability across the good channels
        stat_inds = 1:D.nsamples; 
        stat_inds = stat_inds(1:20:end);
        tmp = std(D(chan_inds(~ch_bad),stat_inds,1),[],1);
        PanWindowData = interp1(linspace(0,1,length(tmp)),tmp,linspace(0,1,D.nsamples));
        PanWindowData(t_bad) = NaN;

        stat_inds = 1:D.nsamples; stat_inds = stat_inds(~t_bad); stat_inds = stat_inds(1:20:end);
        SideWindowData = stat(D(chan_inds,stat_inds,1),[],2);
        SideWindowData = SideWindowData./max(SideWindowData(~ch_bad));
        
    end


    function key_press(~,evnt)
        
        if strcmp(evnt.Key,'rightarrow')
            xspan = get(MainWindow,'XLim');
            set(MainWindow,'XLim',xspan + 0.05*diff(xspan));
            redraw  
        elseif strcmp(evnt.Key,'leftarrow')
            xspan = get(MainWindow,'XLim');
            set(MainWindow,'XLim',xspan - 0.05*diff(xspan));
            redraw      
        end

    end

    function inc_scale(varargin)
        G = G * Ginc;
        redraw
    end

    function dec_scale(varargin)
        G = G / Ginc;
        redraw
    end

    function inc_xwidth(varargin)
        xspan = get(MainWindow,'XLim');
        set(MainWindow,'XLim',mean(xspan) + [-1 1]*diff(xspan))
        redraw
    end

    function dec_xwidth(varargin)
        xspan = get(MainWindow,'XLim');
        set(MainWindow,'XLim',mean(xspan) + 0.25*[-1 1]*diff(xspan))
        redraw
    end

    function mouse_mode(varargin)
        if strcmp(get(uitools.rect,'state'),'on')
            set(hz,'Motion','vertical','Enable','on','RightClickAction','inversezoom')
            setAllowAxesZoom(hz,PanWindow,false);
        else
            set(hz,'Enable','off')
            set(MainFig,'pointer','arrow')
        end
    end

    function ch_bad = get_bad_channels
        % Get bad channels, return them in the currently displayed sort order
        if ~isempty(D.badchannels)
            ch_bad = ismember(chan_labels(chan_inds),D.chanlabels(D.badchannels));
        else
            ch_bad = false(size(chan_inds));
        end
    end

    function t_bad = get_bad_inds(inherited)
        % Return bad time based solely on the BadEpochs variable
        % A kind of stripped-down meeg.badsamples()
        if nargin < 1 || isempty(inherited) 
            inherited = false;
        end

        if inherited
            Bad = InheritedBadEpochs;
        else
            Bad = BadEpochs;
        end
        
        t_bad = false(1,Nsamples);
        for b = 1:numel(Bad.(channel_type))
            if numel(Bad.(channel_type){b}) == 2
                t_bad(D.time>=Bad.(channel_type){b}(1) & D.time<=(Bad.(channel_type){b}(2)+1e-5)) = true;
            end
        end

    end

    function ylims = get_ylims
        ylims = [min(offsets)-2*min(offsets) max(offsets)+2*min(offsets)];
    end

end % OSLview

function [BadEpochs,InheritedBadEpochs] = read_bad_events(D,modality)
    % Set BadEpochs by reading artefact_OSL artefacts from the D object
    % Run once at initialization
    BadEpochs = {};
    InheritedBadEpochs = {};

    ev = D.events;
    if isempty(ev)
        return
    end
    
    % Note - using 'artefact_OSL' here means that we won't interact at all with artefacts identified separately   
    % Otherwise, this function will turn other artefact types into artefact_OSL
    ev = ev(cellfun(@(x) strcmp('artefact_OSL',x),{ev.type}));
    ev2 = ev(ismember({ev.value},{modality})); % These are all the artefact events that apply to the channel types we are inspecting

    for j = 1:numel(ev2)
        BadEpochs{end+1}(1) = ev2(j).time;
        BadEpochs{end}(2) = ev2(j).time + ev2(j).duration - 2/D.fsample;
    end

    if D.montage('getindex')
        m = D.montage('getmontage');
        Dtemp = D.montage('switch',0);
        inherited_chantypes = unique(Dtemp.chantype(find(any(m.tra(D.indchantype(modality),:),1))));
        ev2 = ev(ismember({ev.value},inherited_chantypes)); % These are all the artefact events that apply to the channel types we are inspecting

        for j = 1:numel(ev2)
            InheritedBadEpochs{end+1}(1) = ev2(j).time;
            InheritedBadEpochs{end}(2) = ev2(j).time + ev2(j).duration - 2/D.fsample;
        end
    end

    ev3 = ev(strcmp({ev.value},'all')); % These are all the artefact events that apply to the channel types we are inspecting
    for j = 1:numel(ev3)
        InheritedBadEpochs{end+1}(1) = ev3(j).time;
        InheritedBadEpochs{end}(2) = ev3(j).time + ev3(j).duration - 2/D.fsample;
    end

end

function D = set_bad_events(D,BadEpochs,modality)
    % Save bad epochs using method meeg/events
    BadEvents = struct([]);
    for j = 1:numel(BadEpochs)
        if numel(BadEpochs{j}) == 2
            BadEvents(j).type   = 'artefact_OSL';
            BadEvents(j).value  = modality;
            BadEvents(j).time   = BadEpochs{j}(1);
            BadEvents(j).duration = diff(BadEpochs{j}) + 2/D.fsample; % Need to account for SPM12's weird rounding here
            BadEvents(j).offset = 0;
        end
    end
    
    % Load events
    ev = D.events;
        
    % Remove previous bad epoch events for this
    if isfield(ev,'type')
        to_remove = false(size(ev));
        for j = 1:length(ev)
            if strcmp(ev(j).type,'artefact_OSL') && ismember(ev(j).value,{modality})
                to_remove(j) = 1;
            end
        end
        ev(to_remove) = [];
    end
    
    if ~isempty(BadEvents)
        ev = [ev(:); BadEvents(:)];
    end
    
    % Save new events with previous
    D = events(D,1,ev);
end

function flag = is_multiple_call()
% From File Exchange, Yair Altman
    flag = false;
    % Get the stack
    s = dbstack();
    if numel(s) <= 2
    % Stack too short for a multiple call
    return
    end

    % How many calls to the calling function are in the stack?
    names = {s(:).name};
    TF = strcmp(s(2).name,names);
    count = sum(TF);
    if count>1
    % More than 1
    flag = true;
    end
end
