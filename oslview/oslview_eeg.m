function D = oslview_eeg(D)
% EEG continuous data viewer
% OSLVIEW(D)
% ----------------------------------------------------------------
% D - SPM meeg object
% ----------------------------------------------------------------
% AB 2012 
% Version 1.0 08/08/12 
%
%
% TODO: 
% - check dimensionality
% - zeros marked/detected
% - speed up mag/planar switching
% - go to channel
% - epoched data?
% - filtering/other transformations
% - Allow setting of bad channels by right clicking on the bar chart on the right.
% - Have a reset view button or right click option always available.

disp(['Preparing file ' fullfile(D.path,D.fname)])

% Get directory of the viewer
viewer_dir = strrep(which('oslview'),'oslview.m','');

% Update D with the one on file
datafile = fullfile(D.path,D.fname);
D = spm_eeg_load(datafile);

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
t_bad          = [];
yzoom          = [0 1];
BadEpochs      = {};   
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

% Create plotting windows
MainWindow = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');
PanWindow  = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');
SideWindow = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');

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
uitools.switchchan = uipushtool(uitools.toolbar,    'ClickedCallback',@switch_chantype, 'CData',icons.plan,           'TooltipString','Switch channel type');
uitools.custom     = uipushtool(uitools.toolbar,    'ClickedCallback',@CustomFunction,  'CData',icons.customfunction, 'TooltipString','Apply custom function');

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

drawnow
pointer_wait;


disp('Reading data into memory...')


% Get sample times
Nsamples = D.nsamples;
t = D.time;

% Set downsample factor and initial plotting width
ds = 2;
xwidth = D.fsample*5;
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
channel_type = 'EEG';
% if any(strcmp(D.chantype,'MEGPLANAR'))
%   channel_type = 'MEGPLANAR';
% elseif any(strcmp(D.chantype,'MEGGRAD'))
%   channel_type = 'MEGGRAD';
% end
channel_setup; % will also call redraw & redraw_Sidewindow

% Disable save (and channel switching if CTF)
set(uitools.save,'Enable','off');
if strcmp(channel_type,'MEGGRAD')
  set(uitools.switchchan,'Enable','off');
end

set(MainFig,'visible','on');

pointer_wait;

  function createlayout(varargin) % SETUP FIGURE WINDOWS
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
    
  end


  function redraw % UPDATE PLOTS
    
    ylim_mainwindow = diff(get_ylims).*yzoom + min(get_ylims);
    
    % MainWindow
    cla(MainWindow)
    hold(MainWindow,'on');
    
    % Plot channel signals within current range
    chansig = plot(MainWindow,t(xs),G*ones(size(Nchannels)).*D(chan_inds,xs) + repmat(offsets,size(xs)));
    
    % Add channel label as tag and callback
    for ch = 1:length(chansig)
      set(chansig(ch),'tag',chan_labels{chan_inds(ch)},'ButtonDownFcn',   @line_click);
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
 
    
    % Draw bad events
    redraw_BadEvents(MainWindow,2,'--');
    
    redraw_PanWindow
    
    linkaxes([MainWindow,SideWindow],'y')
    
  end


  function redraw_PanWindow

    PanWindowData_plot = PanWindowData;
    PanWindowData_plot(t_bad) = NaN;
    
    plot(PanWindow,t,PanWindowData_plot)
    axis(PanWindow,'tight')
    ylim = get(PanWindow,'ylim');

    % Add the movable pan box
    box_x = [t(xs(1)) t(xs(end)) t(xs(end)) t(xs(1))];
    box_y = [ylim(1) ylim(1) ylim(2) ylim(2)];
    if ismac
    PanBox = patch(box_x,box_y,'r','parent',PanWindow); %'facealpha',0.5);
    else
    PanBox = patch(box_x,box_y,'r','parent',PanWindow,'facealpha',0.5);
    end
    
    % Draw bad events
    redraw_BadEvents(PanWindow,1,'-');
    
    set(PanWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[]);

  end


  function redraw_SideWindow
    cla(SideWindow)
    hold(SideWindow,'on')
    
    ch_bad = get_bad_channels;
    
    SideWindowData_plot = SideWindowData;
    SideWindowData_plot(ch_bad) = NaN;
    
    % Plot bar data 
    col = colormap(lines); col = col(1:7,:);
    for c = 1:size(col,1)
      chanbar = barh(SideWindow,offsets(c:size(col,1):Nchannels),SideWindowData_plot(c:size(col,1):Nchannels),'facecolor',col(c,:),'barwidth',1/7,'edgecolor','none');  
      % Add context menu to plotted bars
      set(chanbar,'uicontextmenu',ContextMenuS.Menu);
    end
    set(SideWindow,'ylim',get(MainWindow,'ylim'))
    set(SideWindow,'xlim',[min(SideWindowData_plot) max(SideWindowData_plot)])
    set(SideWindow,'xTick',[],'xTicklabel',[],'yTick',[],'yTicklabel',[]);

    redraw_PanWindow
    
  end


  function redraw_BadEvents(ax,lw,ls)
    if ~isempty(BadEpochs)
      for b = 1:numel(BadEpochs)
        line([BadEpochs{b}(1) BadEpochs{b}(1)],get_ylims,'linewidth',lw,'linestyle',ls,'color',[0.1 0.8 0.1],'parent',ax)
        if numel(BadEpochs{b}) == 2
          line([BadEpochs{b}(2) BadEpochs{b}(2)],get_ylims,'linewidth',lw,'linestyle',ls,'color','r','parent',ax)
        end
      end
    end
  end


  function switch_chantype(varargin)
    
    switch channel_type
      case 'MEGPLANAR';
        channel_type = 'MEGMAG';
        set(uitools.switchchan,'CData',icons.mag);
      case 'MEGMAG'
        channel_type = 'MEGPLANAR';
        set(uitools.switchchan,'CData',icons.plan);
    end
    
    channel_setup;
    
  end

  function channel_setup
    pointer_wait;
    % Standard deviation of data for y-axis scaling
    Dsig = std(D(:,1:10:D.nsamples,:),[],2);

    chan_inds = find(strcmp(D.chantype, channel_type));
    
    switch channel_order
      case 'normal'
        chan_inds = sort(chan_inds,'ascend');
      case 'variance'
        [tmp,ix] = sort(Dsig(chan_inds),'ascend');
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
        
    cla(MainWindow), redraw
    redraw_SideWindow
    pointer_wait;
  end


  function ContextMenuOn(src,tmp)
    cp = get(MainWindow,'CurrentPoint');
    
    tp = cp(1,1); % time point
    
    % get channel index from clicked line
    ch = find(strcmp(get(MainWindow,'tag'),chan_labels(chan_inds)));
    
    if isempty(ch) % get closest channel to click location
      ch = cp(1,2);
      [tmp,ch] = min(abs(offsets-ch)); % channel
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

    tmp_line = plot(MainWindow,t(xs),G*ones(size(Nchannels)).*D(chan_inds(ch),xs) + repmat(offsets(ch),size(xs)),'color','k','linewidth',2);
  end


  function ContextMenuOff(varargin)
    if ishandle(tmp_line)
      delete(tmp_line)
    end
  end


 function cb_ContextMenuS(src,tmp)
   
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
   redraw_SideWindow
   pointer_wait;
 end


  function line_click(src,tmp)
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
          [tmp,p1] = min(abs(t-t1(1)));
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
        [tmp,p2] = min(abs(t-t2));
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
        redraw
        calcPlotStats('be')
        redraw_SideWindow
      end
      
            
    elseif strcmp(get(ContextMenuM.MarkEvent,'label'),'Remove Event')  %  Remove existing event   
      get_bad;
      BadEpochs(cellfun(@(x) ~sum(sign(x-t_current)),BadEpochs)) = [];
      set_bad;
      redraw
      calcPlotStats('be')
      redraw_SideWindow
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
    redraw
    redraw_SideWindow
  end


  function calcPlotStats(mode)
    
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
        t_bad = false(1,Nsamples);
        for b = 1:numel(BadEpochs)
          if numel(BadEpochs{b}) == 2
            t_bad(D.time>=BadEpochs{b}(1) & D.time<BadEpochs{b}(2)) = true;
          end
        end
  
        if any(strcmp(mode,{'bc','both'}))
          
          offsets = 3*Dsig;
          offsets(D.badchannels)= nan;
          offsets = offsets(chan_inds);
          offsets(isnan(offsets)) = nanmean(offsets);
          offsets = cumsum(offsets);
          
          
            stat_inds = 1:D.nsamples; stat_inds = stat_inds(1:10:end);
            PanWindowData = std(D(chan_inds(~ch_bad),stat_inds,1));
            PanWindowData = resample(PanWindowData,10,1);
        end
        
        if any(strcmp(mode,{'be','both'}))
            stat_inds = 1:D.nsamples; stat_inds = stat_inds(~t_bad); stat_inds = stat_inds(1:10:end);
            SideWindowData = stat(D(chan_inds,stat_inds,1),[],2);   
            SideWindowData = SideWindowData./max(SideWindowData(~ch_bad));
        end
        
  end


  function key_press(tmp,evnt)
    
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
    [tmp,p1] = min(abs(t-t1));
    
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
    xwidth = xwidth*2;
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


  function mouse_mode(src,tmp)
    if strcmp(get(uitools.rect,'state'),'on')
      set(hz,'Motion','vertical','Enable','on','RightClickAction','inversezoom')
      setAllowAxesZoom(hz,PanWindow,false);
    else
      set(hz,'Enable','off')
      set(MainFig,'pointer','arrow')
    end
  end


  function check_xs
    
    if xwidth > Nsamples
      xwidth = xwidth/2;      
    end
       
    if xs(1) <= 0
      xs = 1:ds:xwidth;
    end
    if xs(end) > Nsamples
      xs = (Nsamples+1:ds:Nsamples+xwidth) - xwidth;
    end
    
   if xwidth < 10
      xwidth = xwidth*2;
      xs = 1:ds:xwidth;
   end
    
  end

   
  function get_bad
    BadEpochs = {};
    Events = D.events;
    for ev = 1:numel(Events)
      if isfield(Events,'type') && strcmp(Events(ev).type,'BadEpoch')
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
        BadEvents(ev).type     = 'BadEpoch';
        BadEvents(ev).value   = ev;
        BadEvents(ev).time     =  BadEpochs{ev}(1);
        BadEvents(ev).duration = diff(BadEpochs{ev});
        BadEvents(ev).offset = 0;
      end
    end
    
    % Load events
    Events = D.events;
        
    % Remove previous bad epoch events
    if isfield(Events,'type')
      Events(strcmp({Events.type},'BadEpoch')) = [];
    end
    
    % Concatenate new and old events
    if size(Events,1) < size(Events,2)
      BadEvents = BadEvents(:);
    end
    if ~isempty(BadEvents)
      Events = [Events(:); BadEvents'];
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


  function save_meeg(varargin)   
    if ~isempty(varargin)
      set(varargin{1},'enable','off');
    end
    
    % Get bad channel and epoch info of current data set:
    get_bad
    ch_bad = get_bad_channels;
    
    D = spm_eeg_load(datafile);
    
    D = badchannels(D,D.meegchannels,0);
    if ~isempty(find(ch_bad,1))
      D = badchannels(D,find(ch_bad),1);
    end
    set_bad;
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
        ch_bad = get_bad_channels;
        
        D = spm_eeg_load(datafile);
        
        set_bad
        D = badchannels(D,D.meegchannels,0);
        if ~isempty(find(ch_bad,1))
          D = badchannels(D,find(ch_bad),1);
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
      copyfile(fullfile(D.path, D.fnamedat),fullfile(Dtmp.path, Dtmp.fnamedat), 'f');
    else
      Dtmp = D;
    end
    
    Dtmp = feval(str2func(strtok(FuncName,'.')),Dtmp);
    
    cd(CurrentPath)
    if strcmp(class(D),'meeg')
      D = Dtmp;
    end
    
    pointer_wait
    channel_setup
    pointer_wait
    
    end
    
    
  end





%   function ToolbarFunc(varargin)
%     SelectedFunction = get(varargin{1},'SelectedItem');  
%     switch SelectedFunction
%       case 'Open custom function...'
%        CurrentPath = pwd;
%        [FuncName,FuncPath] = uigetfile('*.m'); 
%        cd(FuncPath)
%        Data = feval(str2func(strtok(FuncName,'.')),D,Data);
%        cd(CurrentPath)
%        
%        if strcmp(class(Data),'meeg')
%          D = Data;
%          Data = D(:,:,:);
%        end
%        
%       case 'Re-order channels by variance'
%         channel_order = 'variance';
%         channel_setup;
%     end   
%     pointer_wait
%     redraw
%     pointer_wait    
%   end




end % OSLview

