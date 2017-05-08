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

% Update D with the one on file
current_montage = D.montage('getindex');
datafile = fullfile(D.path,D.fname);
D = spm_eeg_load(datafile);
D = montage(D,'switch',current_montage);

% Initialise shared variables
p1=[]; p2=[];
t1=[]; t2=[];
Nchannels      = [];
Dsig           = [];
chan_inds      = [];
chan_labels    = []; 
chancols       = [];
SideWindowData = []; 
PanWindowData  = [];
offsets        = [];
tmp_line       = [];
yzoom          = [0 1];
BadEpochs      = {};  
chanbar = []; % Handles to bars in sidebar
chansig = []; % Handle to line plots in main window
badevents_line = [];
badevents_patch = [];
PanWindow_line = [];
PanWindow_box = [];
get_bad;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create UI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MainFig = figure('Name',                  ['OSLview - ' strtok(D.fname,'.')] ,...
                 'NumberTitle',           'off'     ,...
                 'Menubar',               'none'    ,...
                 'DockControls',          'off'     ,...
                 'WindowButtonDownFcn',   @btn_down ,...
                 'WindowButtonUpFcn',     @btn_up   ,...
                 'KeyPressFcn',           @key_press,...
                 'ResizeFcn',             @resize   ,...
                 'CloseRequestFcn',       @close_fig,...
                 'Visible',               'off');               
           
c = onCleanup(@() delete(MainFig));
             
% Create plotting windows
MainWindow = axes('parent',MainFig, 'units','pixels','XGrid','on');
PanWindow  = axes('parent',MainFig, 'units','pixels');
SideWindow = axes('parent',MainFig, 'units','pixels');

% Load tool icons
icons = load([viewer_dir 'oslview.mat']);

% Create toolbars
uitools.toolbar    = uitoolbar;
uitools.save       = uipushtool(uitools.toolbar,    'ClickedCallback',@save_meeg,       'CData',icons.save,           'TooltipString','Save meeg object');
uitools.expand     = uipushtool(uitools.toolbar,    'ClickedCallback',@inc_xwidth,      'CData',icons.expand,         'TooltipString','Increase time window');    
uitools.shrink     = uipushtool(uitools.toolbar,    'ClickedCallback',@dec_xwidth,      'CData',icons.shrink,         'TooltipString','Decrease time window');
uitools.rect       = uitoggletool(uitools.toolbar,  'ClickedCallback',@mouse_mode,      'CData',icons.rect,           'TooltipString','Zoom to channels');  
uitools.zoomin     = uipushtool(uitools.toolbar,    'ClickedCallback',@inc_scale,       'CData',icons.zoomin,         'TooltipString','Increase scale');           
uitools.zoomout    = uipushtool(uitools.toolbar,    'ClickedCallback',@dec_scale,       'CData',icons.zoomout,        'TooltipString','Decrease scale');  
%uitools.switchchan = uipushtool(uitools.toolbar,    'ClickedCallback',@switch_chantype, 'CData',icons.plan,           'TooltipString','Switch channel type');
uitools.custom     = uipushtool(uitools.toolbar,    'ClickedCallback',@CustomFunction,  'CData',icons.customfunction, 'TooltipString','Apply custom function');

% Create menu for channel selection
uitools.menu_channels = uimenu('label','Channels');
channel_types = unique(D.chantype);
for chtype = 1:length(channel_types)
    uimenu(uitools.menu_channels,'label',channel_types{chtype},'Callback',@switch_chantype);
end
    

% Create Context Menu - Main Window
ContextMenuM.Menu         = uicontextmenu('callback',@ContextMenuOn);
ContextMenuM.ChannelLabel = uimenu(ContextMenuM.Menu, 'label','No Channel Selected','Enable','off');
ContextMenuM.MarkEvent    = uimenu(ContextMenuM.Menu, 'label','Mark Event'        , 'callback',@mark_bad);
ContextMenuM.SetBad       = uimenu(ContextMenuM.Menu, 'label','Set Channel as Bad', 'callback',@bad_channel);
set(MainWindow,'uicontextmenu',ContextMenuM.Menu);

% Create Context Menu - Side Window
ContextMenuS.Menu     = uicontextmenu;
ContextMenuS.Reorder  = uimenu(ContextMenuS.Menu, 'label','Reorder channels by variance','callback',@cb_ContextMenuS);
ContextMenuS.Variance = uimenu(ContextMenuS.Menu, 'label','Variance','Checked','on','callback',@cb_ContextMenuS);
ContextMenuS.Kurtosis = uimenu(ContextMenuS.Menu, 'label','Kurtosis','Checked','off','callback',@cb_ContextMenuS);
set(SideWindow,'uicontextmenu',ContextMenuS.Menu);               

% Setup figure layout and plot data
createlayout

pointer_wait;


% Get sample times
Nsamples = D.nsamples;
t = D.time;

% Set downsample factor and initial plotting width
xwidth = D.fsample*D.time(end);
ds = max([1 fix(xwidth/1000)]);
xs = 1:ds:xwidth;


% Set gain for zooming
G = 1;

% Set zooming factor
Ginc = 2;
 
Layout  = [];
PanBox = [];
ClickedWindow = [];

% Handle for zoom control
hz = zoom(MainFig);

% Set channel type
channel_order = 'normal';
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

pointer_wait;



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
    Layout.MainWindow.Y(1) = 2*Layout.BorderWidth +  Layout.PanWindowWidth;
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
    Layout.SideWindow.Y(1) = 2*Layout.BorderWidth +  Layout.PanWindowWidth;
    Layout.SideWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth;
    Layout.SideWindow.position = [Layout.SideWindow.X(1),Layout.SideWindow.Y(1),abs(diff(Layout.SideWindow.X)),abs(diff(Layout.SideWindow.Y))];

    set(MainWindow,'position',Layout.MainWindow.position);
    set(PanWindow,'position',Layout.PanWindow.position);
    set(SideWindow,'position',Layout.SideWindow.position);
    linkaxes([MainWindow,SideWindow],'y')

  end


  function redraw
    % Update all plots
    
    ylim_mainwindow = diff(get_ylims).*yzoom + min(get_ylims);
    
    % Plot channel signals within current range

    chandata = G*ones(size(Nchannels)).*D(chan_inds,xs) + repmat(offsets(:),size(xs));

    for ch = 1:length(chansig)
      set(chansig(ch),'XData',t(xs),'YData',chandata(ch,:),'LineWidth',0.5,'LineStyle','-','tag',chan_labels{chan_inds(ch)});
    end

    % Set plot limits
    set(MainWindow,'xlim',[t(xs(1)) t(xs(end))])   
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

    for ch = 1:Nchannels
      bar_patch(chanbar(ch),offsets(ch),SideWindowData_plot(ch),range(offsets)/Nchannels,col(ch,:),chan_labels{chan_inds(ch)});
    end
    
    set(SideWindow,'ylim',get(MainWindow,'ylim'))
    set(SideWindow,'xlim',[min(SideWindowData_plot) max(SideWindowData_plot)])
    set(SideWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[]);

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

    pointer_wait;
    % Standard deviation of data
    t_bad = get_bad_inds;
    Dsig_inds = 1:10:D.nsamples; Dsig_inds(t_bad(Dsig_inds)) = [];
    Dsig = std(D(:,Dsig_inds,:),[],2);

    chan_inds = find(strcmp(D.chantype, channel_type));
    
    switch channel_order
      case 'normal'
        chan_inds = sort(chan_inds,'ascend');
      case 'variance'
        [~,ix] = sort(Dsig(chan_inds),'ascend');
        chan_inds = chan_inds(ix);
    end
       
    Nchannels = length(chan_inds);

    % Set channel labels & plot colours
    chan_labels = D.chanlabels;
    chancols = colormap(lines); chancols = chancols(1:7,:);
    chancols = repmat(chancols,ceil(Nchannels/7),1);
    chancols = chancols(1:Nchannels,:);
    
    % Calculate side and pan window statistics
    SideWindowData = [];
    PanWindowData = [];
    calcPlotStats('both')
    
    % Make all of the bars

    cla(MainWindow);
    hold(MainWindow,'on')
    chansig = plot(MainWindow,[0 1],ones(2,Nchannels),'ButtonDownFcn',@line_click);
    badevents_line(1) = plot(MainWindow,[NaN NaN],[NaN NaN],'g','LineWidth',2,'LineStyle','--','HitTest','off','YLimInclude','off');
    badevents_line(2) = plot(MainWindow,[NaN NaN],[NaN NaN],'r','LineWidth',2,'LineStyle','--','HitTest','off','YLimInclude','off');
    badevents_patch(1) = patch(MainWindow,nan(4,1),nan(4,1),'k','LineStyle','none','FaceAlpha',0.1,'HitTest','off','YLimInclude','off');

    cla(SideWindow);
    hold(SideWindow,'on')
    for ch = 1:Nchannels
       chanbar(ch) = patch(SideWindow,zeros(1,4),zeros(1,4),'k','ButtonDownFcn',@line_click,'uicontextmenu',ContextMenuM.Menu);  
    end

    cla(PanWindow);
    hold(PanWindow,'on');
    set(PanWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[]);
    PanWindow_line = plot(PanWindow,NaN,NaN,'b');
    PanWindow_box = patch(NaN,NaN,'r','parent',PanWindow,'facealpha',0.3);
    badevents_line(3) = plot(PanWindow,[NaN NaN],[NaN NaN],'g','LineWidth',1,'LineStyle','-','HitTest','off','YLimInclude','off');
    badevents_line(4) = plot(PanWindow,[NaN NaN],[NaN NaN],'r','LineWidth',1,'LineStyle','-','HitTest','off','YLimInclude','off');
    badevents_patch(2) = patch(PanWindow,nan(4,1),nan(4,1),'k','LineStyle','none','FaceAlpha',0.1,'HitTest','off','YLimInclude','off');

    check_xs; 
    redraw % must check xs still in bounds before redrawing, else may get errors (GC)
    pointer_wait;
  end


  function ContextMenuOn(src,~)
    cp = get(MainWindow,'CurrentPoint');
    
    tp = cp(1,1); % time point
    
    % get channel index from clicked line
    ch = find(strcmp(get(MainWindow,'tag'),chan_labels(chan_inds)));
    
    if isempty(ch) % get closest channel to click location
      ch = cp(1,2);
      [~,ch] = min(abs(offsets-ch)); % channel
    end
    
    set(ContextMenuM.ChannelLabel,'label',chan_labels{chan_inds(ch)});
    
    % Change context menu depending on whether channel is marked as bad
    if isempty(find(D.badchannels==find(strcmp(D.chanlabels,get(ContextMenuM.ChannelLabel,'label'))),1))
      set(ContextMenuM.SetBad,'label','Set Channel as Bad');
    else
      set(ContextMenuM.SetBad,'label','Set Channel as Good');
    end
    
    % Change context menu depending on whether time point with bad epoch
    if any(cellfun(@(x) ~sum(sign(x-tp)),BadEpochs))
      set(ContextMenuM.MarkEvent,'label','Remove Event');
    else
      set(ContextMenuM.MarkEvent,'label','Mark Event'); 
    end

    if isempty(tmp_line)
        tmp_line = plot(MainWindow,t(xs),G*ones(size(Nchannels)).*D(chan_inds(ch),xs) + repmat(offsets(ch),size(xs)),'color','k','linewidth',2);
    else
        delete(tmp_line);
        tmp_line = [];
    end
 end


  function ContextMenuOff(varargin)
    if ~isempty(tmp_line)
        try
            delete(tmp_line)
        end
      tmp_line = [];
    end
  end


   
 function cb_ContextMenuS(src,~)
   
   switch get(src,'label')
     case get(ContextMenuS.Reorder,'label')
       if strcmp(get(ContextMenuS.Reorder,'label'),'Reorder channels by variance');
         channel_order = 'variance';
         set(ContextMenuS.Reorder,'label','Reset channel ordering');
       else
         channel_order = 'normal';
         set(ContextMenuS.Reorder,'label','Reorder channels by variance');
       end
       channel_setup;
     case get(ContextMenuS.Variance,'label')
       set(ContextMenuS.Variance,'Checked','on');
       set(ContextMenuS.Kurtosis,'Checked','off');
     case get(ContextMenuS.Kurtosis,'label')
       set(ContextMenuS.Kurtosis,'Checked','on');
       set(ContextMenuS.Variance,'Checked','off');
   end
   
   pointer_wait;
   calcPlotStats('be');
   redraw
   pointer_wait;
 end


  function line_click(src,~)
     set(MainWindow,'tag',get(src,'tag')); 
  end



  function btn_down(varargin)    
    
    yzoom = (get(MainWindow,'ylim') - min(get_ylims))./(diff(get_ylims));

    set(MainWindow,'tag','');   
    ContextMenuOff % deselect highlighted channel
    
    if strcmp(get(MainFig,'selectionType'),'normal') % if left click
      set(MainFig,'WindowButtonMotionFcn',@btn_move); 
      switch(gca)        
        case MainWindow       
          t1 = get(MainWindow,'CurrentPoint');
          [~,p1] = min(abs(t-t1(1)));
          ClickedWindow = MainWindow;        
        case PanWindow
          ClickedWindow = PanWindow;        
      end
    end
    
  end % btn_down


  function btn_up(varargin)
   
    set(MainFig,'WindowButtonMotionFcn', []);
    if ClickedWindow == PanWindow
      updatePanWindow
    end
    
    ClickedWindow = [];
    
    
  end % btn_up


  function btn_move(varargin)
    
    if ClickedWindow == MainWindow
        t2 = get(MainWindow,'CurrentPoint');
        t2 = t2(1);
        [~,p2] = min(abs(t-t2));
        % Only redraw if change in x position greater than 10% of view
        if abs(t1(1)-t2) > 0.1*diff(get(MainWindow,'xlim'))
          xs = xs + round(p1-p2);
          check_xs
          redraw
        end
    else
      updatePanWindow
    end
    
  end % btn_move


  function mark_bad(varargin)
    
    set(uitools.save,'enable','on');
    t_current = get(MainWindow,'CurrentPoint');
    t_current = t_current(1,1);
    t_window  = get(MainWindow,'xlim');
    
    if strcmp(get(ContextMenuM.MarkEvent,'label'),'Mark Event') %  Create new event
      
      
      % if position of marker is within the first or last 1% of the window
      % then assume it was meant to be at the start or end.
      if t_current < t_window(1) + 0.01*diff(t_window)
        t_current = t_window(1);
      elseif t_current > t_window(1) + 0.95*diff(t_window)
        t_current = t_window(2);
      end

      if isempty(BadEpochs) || length(BadEpochs{end}) == 2 % Create new epoch
        get_bad;
        BadEpochs{end+1}(1) = t_current;
        redraw
      else % End previous epoch
        BadEpochs{end}(2) = t_current;
        set_bad;
        calcPlotStats('both')
        redraw
      end
      
            
    elseif strcmp(get(ContextMenuM.MarkEvent,'label'),'Remove Event')  %  Remove existing event   
      get_bad;
      BadEpochs(cellfun(@(x) ~sum(sign(x-t_current)),BadEpochs)) = [];
      set_bad;
      calcPlotStats('both')
      redraw
    end
  

    
  end


  function bad_channel(varargin)
    set(uitools.save,'enable','on');
    bad_chan = find(strcmp(D.chanlabels,get(ContextMenuM.ChannelLabel,'label')));

    if isempty(find(D.badchannels==bad_chan,1))
      D = badchannels(D,bad_chan,1); 
    else
      D = badchannels(D,bad_chan,0);
    end
        
    calcPlotStats('bc')
    if ~isempty(tmp_line)
      delete(tmp_line);
    end
    redraw
  end


  function calcPlotStats(drawmode)
    
    % Get statistic from side window context menu:
    stat_checked = arrayfun(@(x) get(x,'checked'),get(ContextMenuS.Menu,'children'),'uniformoutput',0);
    stat_labels  = arrayfun(@(x) get(x,'label'),get(ContextMenuS.Menu,'children'),'uniformoutput',0);
    switch(char(stat_labels(strcmp(stat_checked,'on'))))
      case 'Variance', stat = @nanvar;
      case 'Kurtosis', stat = @kurtosis;
    end
    
     
        % Bad Channel indices
        ch_bad = get_bad_channels;
        
        % Bad Sample indices
        t_bad = get_bad_inds;
  
        if any(strcmp(drawmode,{'bc','both'}))
          
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
          
        end
        
        if any(strcmp(drawmode,{'be','both'}))
            stat_inds = 1:D.nsamples; stat_inds = stat_inds(~t_bad); stat_inds = stat_inds(1:20:end);
            SideWindowData = stat(D(chan_inds,stat_inds,1),[],2);   
            SideWindowData = SideWindowData./max(SideWindowData(~ch_bad));
        end
        
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


  function updatePanWindow
    t1 = get(PanWindow,'CurrentPoint');
    t1 = t1(1);
    [~,p1] = min(abs(t-t1));
    
    xs = p1-ceil(xwidth/2)+1:ds:p1+fix(xwidth/2);
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
    xwidth = min(xwidth*2,numel(D.time));
    ds = max([1 fix(xwidth/1000)]);
    xc = round(median(xs));
    xs = xc-ceil(xwidth/2)+1:ds:xc+fix(xwidth/2);
    check_xs
    redraw
  end


  function dec_xwidth(varargin)
    xwidth = fix(xwidth/2);
    ds = max([1 fix(xwidth/1000)]);
    xc = round(median(xs));
    xs = xc-ceil(xwidth/2)+1:ds:xc+fix(xwidth/2);
    check_xs
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


  function check_xs
  % checks validity of xs

  % not all problems are necessarily fixed in one pass. 
  % E.g. if xwidth was several times greater than Nsamples.
  % If we've changed something, let's check again. We should keep tabs on
  % how many times we check
  % GC 2014
  maxNchecks = 15;
  iCheck = 0;
  while iCheck < maxNchecks,
      allChecksOK = run_xs_checks;
      if allChecksOK, 
          return;
      end
      iCheck = iCheck + 1;
  end%while
  % if we're here, have gone through several checks and still not fixed
  xwidth = Nsamples;
  xs = 1:ds:xwidth;

  
      function allChecksOK = run_xs_checks
          allChecksOK = 1;
          
          if xwidth > Nsamples
              xwidth = xwidth/2;
              allChecksOK = 0;
          end
          
          if xs(1) <= 0
              xs = 1:ds:xwidth;
              allChecksOK = 0;
          end
          if xs(end) > Nsamples
              xs = (Nsamples:-ds:Nsamples-xwidth);
              xs = xs(end:-1:1);
              allChecksOK = 0;
          end
          
          if xwidth < 10
              xwidth = xwidth*2;
              xs = 1:ds:xwidth;
              allChecksOK = 0;
          end
      end%run_xs_checks
  end%check_xs

   
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
        BadEvents(ev).type     = 'artefact_OSL';
        BadEvents(ev).value   = 'all';
        BadEvents(ev).time     =  BadEpochs{ev}(1);
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


  function pointer_wait
    if(strcmp(get(MainFig,'pointer'),'watch'))
      set(MainFig,'pointer','arrow');
    else
      set(MainFig,'pointer','watch');
    end
    pause(0)
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
    
    pointer_wait
    channel_setup
    pointer_wait
    
    end
    
    
  end

end % OSLview

