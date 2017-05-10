function D = oslview(D)
	% MEG continuous data viewer
	% OSLVIEW(D)
	% ----------------------------------------------------------------
	% D - SPM meeg object
	% ----------------------------------------------------------------
	% AB 2012
	% Updated 21/10/14
	%
	%
	% TODO:
	% - zeros marked/detected
	% - go to channel
	% - epoched data?
	% - Have a reset view button or right click option always available.
	% - Replace chan_labels with D.chanlabels
	% - Tighter control over SideWindow bars

	% Get directory of the viewer
	viewer_dir = strrep(which('oslview'),'oslview.m','');

	% Initialise shared variables
	drag = struct('state',0); % Struct to track the state of the main window drag
	chan_time = []; % Array of time points for plots on main window
	chan_data = []; % Matrix of value points for plots on main window
	chan_inds = []; % Array specifying order in which to display the components
	chan_labels = [];

	selection = struct('line',[],'bar',[]); % Temporary lines for context menu

	Nchannels = [];
	Dsig = [];
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
	get_bad;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create UI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	MainFig = figure('Name',['OSLview - ' strtok(D.fname,'.')],...
						'NumberTitle','off',...
						'Menubar','none',...
						'DockControls','off',...
						'WindowButtonDownFcn',@btn_down	,...
						'WindowButtonUpFcn',@btn_up	,...
						'KeyPressFcn',@key_press,...
						'ResizeFcn',@resize	,...
						'CloseRequestFcn',@close_fig,...
                        'Units','pixels',...
						'Visible','off');
			
	c = onCleanup(@() delete(MainFig));
				
	% Create plotting windows
	MainWindow = axes('parent',MainFig, 'units','pixels','XGrid','on');
	PanWindow	= axes('parent',MainFig, 'units','pixels');
	SideWindow = axes('parent',MainFig, 'units','pixels');

	% Load tool icons
	icons = load([viewer_dir 'oslview.mat']);

	% Create toolbars
	uitools.toolbar		= uitoolbar;
	uitools.save				= uipushtool(uitools.toolbar,		'ClickedCallback',@save_meeg,				'CData',icons.save,						'TooltipString','Save meeg object');
	uitools.expand			= uipushtool(uitools.toolbar,		'ClickedCallback',@inc_xwidth,			'CData',icons.expand,					'TooltipString','Increase time window');
	uitools.shrink			= uipushtool(uitools.toolbar,		'ClickedCallback',@dec_xwidth,			'CData',icons.shrink,					'TooltipString','Decrease time window');
	uitools.rect				= uitoggletool(uitools.toolbar,	'ClickedCallback',@mouse_mode,			'CData',icons.rect,						'TooltipString','Zoom to channels');
	uitools.zoomin			= uipushtool(uitools.toolbar,		'ClickedCallback',@inc_scale,				'CData',icons.zoomin,					'TooltipString','Increase scale');
	uitools.zoomout		= uipushtool(uitools.toolbar,		'ClickedCallback',@dec_scale,				'CData',icons.zoomout,				'TooltipString','Decrease scale');
	%uitools.switchchan = uipushtool(uitools.toolbar,		'ClickedCallback',@switch_chantype, 'CData',icons.plan,						'TooltipString','Switch channel type');
	uitools.custom			= uipushtool(uitools.toolbar,		'ClickedCallback',@CustomFunction,	'CData',icons.customfunction, 'TooltipString','Apply custom function');

	% Create menu for channel selection
	uitools.menu_channels = uimenu('label','Channels');
	channel_types = unique(D.chantype);
	for chtype = 1:length(channel_types)
		uimenu(uitools.menu_channels,'label',channel_types{chtype},'Callback',@switch_chantype);
	end
		

	% Create Context Menu - Main Window
	ContextMenuM.Menu					= uicontextmenu('callback',@ContextMenuOn);
	ContextMenuM.ChannelLabel = uimenu(ContextMenuM.Menu, 'label','No Channel Selected','Enable','off');
	ContextMenuM.MarkEvent		= uimenu(ContextMenuM.Menu, 'label','Mark Event'				, 'callback',@mark_bad);
	ContextMenuM.SetBad				= uimenu(ContextMenuM.Menu, 'label','Set Channel as Bad', 'callback',@bad_channel);
	set(MainWindow,'uicontextmenu',ContextMenuM.Menu);

	% Create Context Menu - Side Window axis
	ContextMenuS.Menu			= uicontextmenu;
	ContextMenuS.Reorder	= uimenu(ContextMenuS.Menu, 'label','Sort by metric','callback',@cb_ContextMenuS);
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

	Layout	= [];
	PanBox = [];
	ClickedWindow = [];

	% Handle for zoom control
	hz = zoom(MainFig);

	% Set channel type
	channel_order = 'normal'; % Default is to use normal ordering (by index)
	if any(strcmp(D.chantype,'MEGMAG')) && any(strcmp(D.chantype,'MEGPLANAR')) % ELEKTA
		channel_type = 'MEGPLANAR';
	elseif any(strcmp(D.chantype,'MEGMAG')) % 4D
		channel_type = 'MEGMAG';
	elseif any(strcmp(D.chantype,'MEGGRAD')) % CTF
		channel_type = 'MEGGRAD';
	elseif any(strcmp(D.chantype,'MEG'))
		channel_type = 'MEG';
	else % catch unrecognised case
		channel_type = D.chantype(1);
	end
	channel_setup; % will also call redraw & redraw_Sidewindow

	% Disable save (and channel switching if CTF)
	set(uitools.save,'Enable','off');

	set(MainFig,'visible','on');

	drawnow
	warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
	jFig = get(handle(MainFig),'JavaFrame');
	jFig.setMaximized(true);
	warning on MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame

	%drag_listener1 = addlistener(MainWindow,'XLim','PostSet',@redraw);
	set(PanWindow,'ButtonDownFcn',@pan_jump)


	h_pan = pan(MainFig);
	setAllowAxesPan(h_pan,PanWindow,false)
	setAllowAxesPan(h_pan,SideWindow,false)
	setAllowAxesPan(h_pan,MainWindow,true)
	setAxesPanConstraint(h_pan,MainWindow,'x')
	h_pan.Enable = 'off';
	% h_pan.ButtonDownFilter = @pan_filter
	h_pan.ActionPreCallback = @plot_low_res;
	h_pan.ActionPostCallback = @redraw;



	function flag = pan_filter(obj,event)
		if obj == PanWindow || obj == SideWindow
		flag = true;
		elseif ~strcmp(get(MainFig,'selectionType'),'normal'); % Pass through all right clicks
		flag = true;
		else
		flag = false;
		end
	end



	% drag_listener2 = addlistener(MainWindow,'XLim','PreSet',@validate_xlim);

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
		Layout.MainWindow.Y(1) = 2*Layout.BorderWidth +	Layout.PanWindowWidth;
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
		Layout.SideWindow.Y(1) = 2*Layout.BorderWidth +	Layout.PanWindowWidth;
		Layout.SideWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth;
		Layout.SideWindow.position = [Layout.SideWindow.X(1),Layout.SideWindow.Y(1),abs(diff(Layout.SideWindow.X)),abs(diff(Layout.SideWindow.Y))];

		set(MainWindow,'position',Layout.MainWindow.position,'XLim',D.time([1 end]));
		set(PanWindow,'position',Layout.PanWindow.position);
		set(SideWindow,'position',Layout.SideWindow.position);
		box(MainWindow,'on')
		box(PanWindow,'on')
		box(SideWindow,'on')

	end

	function pan_jump(ax,hit)
		xspan = get(MainWindow,'XLim')
		set(MainWindow,'XLim',hit.IntersectionPoint(1)+[-0.5 0.5]*diff(xspan));
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


		% set(MainWindow,'xlim',[t(xs(1)) t(xs(end))])
		set(MainWindow,'ylim',ylim_mainwindow)
		set(MainWindow,'ytick',[],'yticklabel',[])

		% Add context menu to plotted signals and label with channel name
		set(chansig,'uicontextmenu',ContextMenuM.Menu);
		
		% Replace bad channels with dashed lines
		ch_bad = get_bad_channels;
		
		set(chansig(ch_bad),'linestyle','--','linewidth',1);

		%% PAN WINDOW
		PanWindowData_plot = PanWindowData;
		PanWindowData_plot(get_bad_inds) = NaN;
		
		set(PanWindow_line,'XData',t,'YData',PanWindowData_plot);
		axis(PanWindow,'tight')
		ylim = get(PanWindow,'ylim');

		box_x = [t(xs(1)) t(xs(end)) t(xs(end)) t(xs(1))];
		box_y = [ylim(1) ylim(1) ylim(2) ylim(2)];
		set(PanWindow_box,'XData',box_x,'YData',box_y);

		%% SIDE WINDOW
		
		% Read offsets for the bar positions from elsewhere as well
		ch_bad = get_bad_channels;
		SideWindowData_plot = SideWindowData;
		SideWindowData_plot(ch_bad) = NaN;
		
		% Plot bar data
		col = colormap(lines);
		col = col(1:7,:);
		col = repmat(col,ceil(Nchannels/7),1);

		% Because offsets are set to cumsum when they are computed, they are guaranteed to be sorted
		% So they can be replaced here with just simple offsets
		sidebar_offsets = (1:length(offsets))-0.5;
		for ch = 1:Nchannels
			bar_patch(chanbar(ch),sidebar_offsets(ch),SideWindowData_plot(ch),1,col(ch,:),chan_labels{chan_inds(ch)});
		end
		
		set(SideWindow,'ylim',[0 length(offsets)])
		set(SideWindow,'xlim',[min(SideWindowData_plot) max(SideWindowData_plot)])
		set(SideWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[]);

		if strcmp(get(ContextMenuS.Variance,'Checked'),'on')
			xlabel(SideWindow,'Variance');
		else
			xlabel(SideWindow,'Kurtosis');
		end

		%% BAD EVENTS
		yl = get_ylims;
		green_x = NaN;
		green_y = NaN;
		red_x = NaN;
		red_y = NaN;
		patch_x = [];
		patch_y = [];

		for j = 1:length(BadEpochs)
			green_x = [green_x BadEpochs{j}(1) BadEpochs{j}(1) NaN];
			green_y = [green_y yl NaN];
			if numel(BadEpochs{j}) == 2
			red_x = [red_x BadEpochs{j}(2) BadEpochs{j}(2) NaN];
			red_y = [red_y yl NaN];
			patch_x(:,end+1) = [BadEpochs{j}(1);BadEpochs{j}(2);BadEpochs{j}(2);BadEpochs{j}(1)];
			patch_y(:,end+1) = [yl(1);yl(1);yl(2);yl(2)];
			end
		end

		set(badevents_line([1,3]),'XData',green_x,'YData',green_y);
		set(badevents_line([2,4]),'XData',red_x,'YData',red_y);
		set(badevents_patch,'XData',patch_x,'YData',patch_y);
		
	end

	function bar_patch(h_patch,offset,value,width,c,tag)
		% Helper function to use a patch object as a fast horizontal bar
		xv = [0 value];
		yv = offset + [-1 1]*width/2;
		set(h_patch,'YData',[yv(1) yv(1) yv(2) yv(2)],'XData',[xv(1) xv(2) xv(2) xv(1)],'FaceColor',c,'LineStyle','none','tag',tag)
	end

	function switch_chantype(src,~)
		
		channel_type = get(src,'label');
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
		chan_labels = D.chanlabels;
		chancols = colormap(lines); 
		chancols = chancols(1:7,:);
		chancols = repmat(chancols,ceil(Nchannels/7),1);
		chancols = chancols(1:Nchannels,:);
		
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
		badevents_patch(1) = patch(MainWindow,nan(4,1),nan(4,1),'k','LineStyle','none','FaceAlpha',0.1,'HitTest','off','YLimInclude','off');

		cla(SideWindow);
		hold(SideWindow,'on')
		for ch = 1:Nchannels
			chanbar(ch) = patch(SideWindow,zeros(1,4),zeros(1,4),'k','uicontextmenu',ContextMenuSP.Menu);
		end

		cla(PanWindow);
		hold(PanWindow,'on');
		set(PanWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[],'ButtonDownFcn',@pan_jump);
		PanWindow_line = plot(PanWindow,NaN,NaN,'b','HitTest','off');
		PanWindow_box = patch(NaN,NaN,'r','parent',PanWindow,'facealpha',0.3,'HitTest','off');
		badevents_line(3) = plot(PanWindow,[NaN NaN],[NaN NaN],'g','LineWidth',1,'LineStyle','-','HitTest','off','YLimInclude','off');
		badevents_line(4) = plot(PanWindow,[NaN NaN],[NaN NaN],'r','LineWidth',1,'LineStyle','-','HitTest','off','YLimInclude','off');
		badevents_patch(2) = patch(PanWindow,nan(4,1),nan(4,1),'k','LineStyle','none','FaceAlpha',0.1,'HitTest','off','YLimInclude','off');

		redraw

	end


	function ContextMenuOn(src,~)
		
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
		if any(cellfun(@(x) ~sum(sign(x-cp(1))),BadEpochs))
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
		set(MainWindow,'tag','');
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

		drag = struct;
		drag.initial_point = get(MainWindow,'CurrentPoint');
		drag.initial_xlim = get(MainWindow,'XLim');
		drag.state = 2; % 0 - no drag, 1 - drag in progress, 2 - drag armed

		set(MainWindow,'tag','');
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
		set(MainFig,'WindowButtonMotionFcn',[]);
		if drag.state == 1
			redraw
		end
		drag.state = 0;
	end

	function mark_bad(varargin)
		
		set(uitools.save,'enable','on');
		t_current = get(MainWindow,'CurrentPoint');
		t_current = t_current(1,1);
		t_window	= get(MainWindow,'xlim');
		
		if strcmp(get(ContextMenuM.MarkEvent,'label'),'Mark Event') %	Create new event
			% if position of marker is within the first or last 1% of the window
			% then assume it was meant to be at the start or end.
			if t_current < t_window(1) + 0.01*diff(t_window)
				t_current = t_window(1);
			elseif t_current > t_window(1) + 0.99*diff(t_window)
				t_current = t_window(2);
			end

			if isempty(BadEpochs) || length(BadEpochs{end}) == 2 % Create new epoch
				get_bad;
				BadEpochs{end+1}(1) = t_current;
				redraw
			else % End previous epoch
				BadEpochs{end}(2) = t_current;
				set_bad;
				calcPlotStats
				redraw
			end		
		elseif strcmp(get(ContextMenuM.MarkEvent,'label'),'Remove Event')	%	Remove existing event
			get_bad;
			BadEpochs(cellfun(@(x) ~sum(sign(x-t_current)),BadEpochs)) = [];
			set_bad;
			calcPlotStats
			redraw
		end
	end


	function bad_channel(varargin)
		% Read the bad channel from the parent 
		label = get(findobj(get(varargin{1},'Parent'),'Enable','off'),'Label'); % Get the parent ContextMenu, find the entry that's greyed out, and return its label

		set(uitools.save,'enable','on');
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
		t_bad = get_bad_inds;

		
		Dsig_inds = 1:20:D.nsamples; Dsig_inds(t_bad(Dsig_inds)) = [];
		Dsig = std(D(:,Dsig_inds,:),[],2);
		offsets = 3*iqr(D(:,Dsig_inds,:),2);
		offsets(offsets==0) = eps;
		offsets(D.badchannels) = nan;
		offsets = offsets(chan_inds);
		offsets(isnan(offsets)) = nanmean(offsets);
		offsets = cumsum(offsets);
	
	
		stat_inds = 1:D.nsamples; stat_inds = stat_inds(1:20:end);
		tmp = std(D(chan_inds(~ch_bad),stat_inds,1));
		PanWindowData = interp1(linspace(0,1,length(tmp)),tmp,linspace(0,1,D.nsamples));
		

		stat_inds = 1:D.nsamples; stat_inds = stat_inds(~t_bad); stat_inds = stat_inds(1:20:end);
		SideWindowData = stat(D(chan_inds,stat_inds,1),[],2);
		SideWindowData = SideWindowData./max(SideWindowData(~ch_bad));
		
	end


	function key_press(~,evnt)
		
		if strcmp(evnt.Key,'rightarrow') || strcmp(evnt.Key,'leftarrow')
			key_scroll(evnt.Key)
		end

	end % key_press


	function key_scroll(key)
		if strcmp(key,'rightarrow')
			scroll_dir = 1;
		else
			scroll_dir = -1;
		end
		
		xs = xs + scroll_dir*round(0.1*diff(xs([1 end])));
		check_xs
		redraw
	end

	function resize(varargin)
		createlayout
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


	function mouse_mode(src,~)
		if strcmp(get(uitools.rect,'state'),'on')
			set(hz,'Motion','vertical','Enable','on','RightClickAction','inversezoom')
			setAllowAxesZoom(hz,PanWindow,false);
		else
			set(hz,'Enable','off')
			set(MainFig,'pointer','arrow')
		end
	end

	function get_bad
		BadEpochs = {};
		Events = D.events;
		for ev = 1:numel(Events)
			if isfield(Events,'type') && strcmp(Events(ev).type,'artefact_OSL')
				BadEpochs{end+1}(1) = Events(ev).time;
				BadEpochs{end}(2) = Events(ev).time + Events(ev).duration;
			end
		end
	end

	function set_bad
		% Save bad epochs using method meeg/events
		BadEvents = struct([]);
		for ev = 1:numel(BadEpochs)
			if numel(BadEpochs{ev} == 2)
			BadEvents(ev).type			= 'artefact_OSL';
			BadEvents(ev).value		= 'all';
			BadEvents(ev).time			=	BadEpochs{ev}(1);
			BadEvents(ev).duration = diff(BadEpochs{ev});
			BadEvents(ev).offset = 0;
			end
		end
		
		% Load events
		Events = D.events;
			
		% Remove previous bad epoch events
		if isfield(Events,'type')
			Events(strcmp({Events.type},'artefact_OSL')) = [];
		end
		
		% Concatenate new and old events
		if size(Events,1) < size(Events,2)
			BadEvents = BadEvents(:);
		end
		if ~isempty(BadEvents)
			Events = [Events(:); BadEvents(:)];
		end
		
		% Save new events with previous
		D = events(D,1,Events);
	end

	function ch_bad = get_bad_channels
		if ~isempty(D.badchannels)
			ch_bad = ismember(chan_labels(chan_inds),D.chanlabels(D.badchannels));
		else
			ch_bad = false(size(chan_inds));
		end
	end

	function t_bad = get_bad_inds
		t_bad = false(1,Nsamples);
		for b = 1:numel(BadEpochs)
			if numel(BadEpochs{b}) == 2
				t_bad(D.time>=BadEpochs{b}(1) & D.time<=BadEpochs{b}(2)) = true;
			end
		end
	end

	function save_meeg(varargin)
		if ~isempty(varargin)
			set(varargin{1},'enable','off');
		end
		
		% Get bad channel and epoch info of current data set:
		get_bad
		ch_bad = D.badchannels;
		
		D = spm_eeg_load(datafile);
		
		D = badchannels(D,D.indchantype('MEEG'),0);
		if ~isempty(ch_bad)
			D = badchannels(D,ch_bad,1);
		end
		set_bad;
		D = montage(D,'switch',current_montage);

		save(D);
	end



	function close_fig(varargin)
		if strcmp(get(uitools.save,'enable'),'on')
			save_flag = questdlg('Save MEEG figure before closing?');
			switch save_flag
			case 'Yes'
				save_meeg;
				delete(MainFig)
			case {'Cancel',''}
				return
			case 'No'
				delete(MainFig)
			end
		else
			delete(MainFig)

		end
	end



	function ylims = get_ylims
		ylims = [min(offsets)-2*min(offsets) max(offsets)+2*min(offsets)];
	end

	function CustomFunction(varargin)
		
		try
			
			if ~strcmp(datafile,fullfile(D.path,D.fname))
				choice = questdlg('Use current or original data?','OSLview','Current','Original','Cancel','Original');
				switch choice
					case 'Original'
						% Get bad channel and epoch info of current data set:
						get_bad
						ch_bad = D.badchannels;
						D = spm_eeg_load(datafile);
				
						set_bad
						D = badchannels(D,D.indchantype('MEEG'),0);
						if ~isempty(ch_bad)
							D = badchannels(D,ch_bad,1);
						end
						save(D);
					case 'Cancel'
						return
				end
			end
			
			CurrentPath = fileparts(mfilename('fullpath'));
			[FuncName,FuncPath] = uigetfile([CurrentPath,'/*.m']);
			cd(FuncPath)
			
			if ~strcmp(fullfile(D.path,D.fname),[CurrentPath,'/Dtmp.mat'])
				Dtmp = clone(D,[CurrentPath,'/Dtmp']);
				copyfile(D.fnamedat,Dtmp.fnamedat, 'f');
			else
				Dtmp = D;
			end
			
			Dtmp = feval(str2func(strtok(FuncName,'.')),Dtmp);
			
			cd(CurrentPath)
			if isa(D,'meeg')
				D = Dtmp;
			end
			
			channel_setup
		
		end
		
		
	end

end % OSLview

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
