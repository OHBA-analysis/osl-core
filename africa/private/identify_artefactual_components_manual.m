function bad_components = identify_artefactual_components_manual(D,tc,topos,metrics,bad_components)
% Note - tc is passed in with NaNs already in place
% But the artefact channels need to have NaNs added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create UI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MainFig = figure('Name',            ['AFRICA: ' D.fname], ...
                 'NumberTitle',     'off'               , ...
                 'Menubar',         'none'              , ...
                 'DockControls',    'on'                , ...
                 'KeyPressFcn',     @key_press          , ...
                 'Interruptible',   'off'               , ... 
                 'BusyAction',      'queue'            , ...
                 'Visible',         'off'               , ...
                 'ResizeFcn',       @createlayout); % ,...

goodcolor    = [ 48   128    20] / 255;
badcolor     = [204     0     0] / 255;
currentcolor = [100   149   237] / 255;
covcolor     = [121   121   121] / 255;
FONTSIZE     = 14;

metric_names = fieldnames(metrics);
current_metric = 1;

[~,sorted_comps] = sort(metrics.(metric_names{current_metric}).value,'descend');

current_comp_ind = 1;
current_comp     = sorted_comps(current_comp_ind);


% Create plotting windows
tICWindow    = axes('parent',MainFig, 'units','pixels');
specWindow   = axes('parent',MainFig, 'units','pixels');
metricWindow = axes('parent',MainFig, 'units','pixels');
covWindow    = axes('parent',MainFig, 'units', 'pixels');
topoWindow = [];

% If topos not precomputed, then compute them now
if isempty(topos)
    fprintf('No precomputed topos present, computing them now...\n')
    topos = [];
    modalities = D.ica.modalities; 
    for m = 1:numel(modalities)
        topos = [topos component_topoplot(D,D.ica.sm,modalities(m))];
    end
else
    topos = D.ica.topos;
end


for m = 1:size(topos,2)
    topoWindow(m) = axes('parent',MainFig, 'units','pixels');
end

% icons
[icon_good,icon_bad,icon_zoom] = load_icons;

uitools.toolbar    = uitoolbar;
uitools.setgood    = uipushtool(uitools.toolbar,    'ClickedCallback',@setgood, 'CData',icon_good, 'TooltipString','Set component as good');
uitools.setbad     = uipushtool(uitools.toolbar,    'ClickedCallback',@setbad,  'CData',icon_bad,  'TooltipString','Set component as bad');    
uitools.zoom       = uitoggletool(uitools.toolbar,  'ClickedCallback',@cb_zoom, 'CData',icon_zoom,  'TooltipString','Set component as bad');    
uitools.metrics = uicontrol('Style', 'popup', 'String', metric_names, 'Position', [1 1 1 1], 'Callback', @switchmetric);   

% Create context menu for side window
metricContext.menu   = uicontextmenu;
metricContext.switch = uimenu(metricContext.menu, 'label','Reorder channels on metric switch','Checked','on','callback',@cb_metricContext);
set(metricWindow,'uicontextmenu',metricContext.menu);      

drawnow
redraw

try
    warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
    jFig = get(handle(MainFig),'JavaFrame');
    jFig.setMaximized(true);
    warning on MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
catch
    fprintf(2,'Could not maximize window - code may need to be updated\n')
end

set(MainFig,'Visible','on')
uiwait(MainFig)

    function createlayout(varargin) % SETUP FIGURE WINDOWS
        Layout.MainFig.position = get(MainFig,'position');
        Layout.FigWidth         = Layout.MainFig.position(3);
        Layout.FigHeight        = Layout.MainFig.position(4);
        Layout.BorderWidth      = 20;
        Layout.LabelWidth       = 20;
        Layout.MenuHeight       = 20;
        
        vAxisSpace                = Layout.FigHeight - 2*Layout.BorderWidth - 8*Layout.LabelWidth;
        
        Layout.tICWindowHeight    = fix(0.35*vAxisSpace);
        Layout.covWindowHeight    = fix(0.25*vAxisSpace);
        
        Layout.specWindowHeight   = fix(0.4*vAxisSpace);
        Layout.specWindowWidth    = fix(0.45*(Layout.FigWidth  - 3*Layout.BorderWidth - 3*Layout.LabelWidth));
        
        Layout.topoWindowHeight   = fix(0.4*vAxisSpace);
        Layout.topoWindowWidth    = fix(0.45*(Layout.FigWidth  - 3*Layout.BorderWidth - 3*Layout.LabelWidth));
        
        Layout.metricWindowHeight = Layout.FigHeight           - 2*Layout.BorderWidth - 2*Layout.LabelWidth;
        Layout.metricWindowWidth  = fix(0.1*(Layout.FigWidth   - 3*Layout.BorderWidth - 3*Layout.LabelWidth));

        
        
        % tIC Window
        Layout.tICWindow.X(1) = Layout.BorderWidth + 2*Layout.LabelWidth;
        Layout.tICWindow.X(2) = Layout.FigWidth - Layout.BorderWidth - Layout.LabelWidth - Layout.metricWindowWidth;
        Layout.tICWindow.Y(1) = Layout.FigHeight - Layout.BorderWidth - Layout.LabelWidth - Layout.tICWindowHeight;
        Layout.tICWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth - Layout.LabelWidth;
        Layout.tICWindow.position = [Layout.tICWindow.X(1),Layout.tICWindow.Y(1) ,abs(diff(Layout.tICWindow.X)),abs(diff(Layout.tICWindow.Y))];
        
        
        % covariate Window
        Layout.covWindow.X(1) = Layout.BorderWidth + 2*Layout.LabelWidth;
        Layout.covWindow.X(2) = Layout.FigWidth - Layout.BorderWidth - Layout.LabelWidth - Layout.metricWindowWidth;
        Layout.covWindow.Y(1) = Layout.BorderWidth + 3*Layout.LabelWidth + Layout.specWindowHeight;
        Layout.covWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth - 4*Layout.LabelWidth - Layout.tICWindowHeight;
        Layout.covWindow.position = [Layout.covWindow.X(1),Layout.covWindow.Y(1) ,abs(diff(Layout.covWindow.X)),abs(diff(Layout.covWindow.Y))];
        
        
        % spectrum window
        Layout.specWindow.X(1) = Layout.BorderWidth + 2*Layout.LabelWidth;
        Layout.specWindow.X(2) = Layout.specWindow.X(1) + Layout.specWindowWidth;
        Layout.specWindow.Y(1) = Layout.BorderWidth + Layout.LabelWidth;
        Layout.specWindow.Y(2) = Layout.BorderWidth + Layout.LabelWidth + Layout.specWindowHeight;
        Layout.specWindow.position = [Layout.specWindow.X(1),Layout.specWindow.Y(1) ,abs(diff(Layout.specWindow.X)),abs(diff(Layout.specWindow.Y))];
        
        
        % topoplot window
        Layout.topoWindow.X(1) = Layout.BorderWidth + 3*Layout.LabelWidth + Layout.specWindowWidth;
        Layout.topoWindow.X(2) = Layout.FigWidth - Layout.BorderWidth - Layout.LabelWidth - Layout.metricWindowWidth;

        Layout.topoWindow.Y(1) = Layout.BorderWidth + Layout.LabelWidth;
        Layout.topoWindow.Y(2) = Layout.topoWindow.Y(1) + Layout.topoWindowHeight;
        Layout.topoWindow.position = [Layout.topoWindow.X(1),Layout.topoWindow.Y(1),abs(diff(Layout.topoWindow.X)),abs(diff(Layout.topoWindow.Y))];
        
        
        % metric Window
        Layout.metricWindow.X(1) = Layout.BorderWidth + 4*Layout.LabelWidth + Layout.specWindowWidth + Layout.topoWindowWidth;
        Layout.metricWindow.X(2) = Layout.FigWidth - Layout.BorderWidth;
        Layout.metricWindow.Y(1) = Layout.BorderWidth + Layout.LabelWidth;
        Layout.metricWindow.Y(2) = Layout.FigHeight - Layout.BorderWidth - Layout.LabelWidth;
        Layout.metricWindow.position = [Layout.metricWindow.X(1),Layout.metricWindow.Y(1) ,abs(diff(Layout.metricWindow.X)),abs(diff(Layout.metricWindow.Y))];
        
        


        if all([Layout.tICWindow.position([3,4]) Layout.covWindow.position([3,4]) Layout.specWindow.position([3,4]) Layout.topoWindow.position([3,4])] > 0)
            set(tICWindow,   'position',Layout.tICWindow.position);
            set(specWindow,  'position',Layout.specWindow.position);            
            set(metricWindow,'position',Layout.metricWindow.position);
            set(covWindow,   'position',Layout.covWindow.position);
            if length(topoWindow) == 2
                toposub1 = fix(Layout.topoWindow.position .* [1,1,1/length(topoWindow),1]);
                toposub2 = [toposub1(1) + toposub1(3),toposub1(2:4)];
                set(topoWindow(1),'position',toposub1);
                set(topoWindow(2),'position',toposub2);
            else
                set(topoWindow,'position',Layout.topoWindow.position);
            end
        end
        
        % metrics dropdown menu layout
        set(uitools.metrics,'position',[Layout.metricWindow.position(1) Layout.metricWindow.Y(2) Layout.metricWindow.position(3) Layout.MenuHeight])
        linkaxes([tICWindow,covWindow],'x');

    end




    function redraw % UPDATE PLOTS
          
        drawnow expose
        
        t = (1:size(tc,2))./D.fsample;
        
        if ismember(current_comp,bad_components)
            thistcplotcolor = badcolor;
        else
            thistcplotcolor = goodcolor; %tcplotcolor;
        end
             
        
        % Redraw topo window
        for m = 1:length(topoWindow)

            % clear axis contents
            cla(topoWindow(m));            
            struct2handle(topos(current_comp,m).children,topoWindow(m)); 
            
            % Set colormap limits
            contourGroupInd = strcmpi('specgraph.contourgroup', {topos(current_comp,m).children(:).type});
            cmax = max(abs(topos(current_comp,m).children(contourGroupInd).properties.ZData(:)));
            
            struct2handle(topos(current_comp,m).children,topoWindow(m)); 
            axis(topoWindow(m),'image')
            set(topoWindow(m),'Visible','off')
            set(topoWindow(m), 'CLim', cmax*[-1 1]);
            set(MainFig,'CurrentAxes',topoWindow(m)); colormap(bluewhitered);

        end
        
        % Redraw spectrum window
        axes(specWindow); 
        
        % plot(component_f,component_P(:,current_comp));
        % set(specWindow,'YScale','log','XLim',[0 D.fsample/2])
        % ylabel('Power spectral density')
        % xlabel('Frequency (Hz)');
        
        pwelch(tc(current_comp,~isnan(tc(current_comp,:))),1024,512,1024,D.fsample);
        
        set(findobj(specWindow,'type','line'),'color', thistcplotcolor,'linewidth',2);
        title(specWindow,'')
        tidyAxes(specWindow, FONTSIZE);
        
        % redraw covariate window
        axes(covWindow);
        if isfield(metrics.(metric_names{current_metric}), 'chanind')
            cov_tc = D(metrics.(metric_names{current_metric}).chanind,:,:);
            cov_tc = reshape(cov_tc,1,D.nsamples*D.ntrials);
            cov_tc(event_to_sample(D,'artefact_OSL',D.chantype(metrics.(metric_names{current_metric}).chanind))) = NaN;
            plot(t(:), cov_tc,'color',covcolor);
        else
            plot(t, zeros(size(t)), 'color', covcolor);
        end%if
        
        axis tight
        set(covWindow,'xlim',[t(1) t(end)]);
        set(covWindow,'ylim',max(abs(get(covWindow,'ylim')))*[-1 1]);
        xlabel('Time (s)', 'FontSize', FONTSIZE);
        tidyAxes(covWindow, FONTSIZE)
        
        % Redraw tIC window
        axes(tICWindow); 
        cla(tICWindow);
        plot(t,tc(current_comp,:),'color',thistcplotcolor), 
        set(tICWindow,'ylim',max(abs(get(tICWindow,'ylim')))*[-1 1]);
        set(tICWindow,'xlim',[t(1) t(end)]);
        tidyAxes(tICWindow, FONTSIZE);

        
        % Redraw metric window
        axes(metricWindow), cla(metricWindow), hold(metricWindow,'on');
        
        barMetric = metrics.(metric_names{current_metric}).value(sorted_comps); 
        barInd = 1:length(barMetric);
        
        [goodbars,badbars,currentbar] = deal(barMetric);
        
        goodbars(   ismember(sorted_comps,bad_components)) = nan;
        badbars(   ~ismember(sorted_comps,bad_components)) = nan;
        currentbar(~ismember(sorted_comps,current_comp))   = nan;

        h_goodbars      = barh(barInd,goodbars);
        h_badbars       = barh(barInd,badbars);
        h_currentbar    = barh(barInd,currentbar);
        
        set(h_goodbars,  'FaceColor', goodcolor,    'EdgeColor', 'none','HitTest','off');
        set(h_badbars,   'FaceColor', badcolor,     'EdgeColor', 'none','HitTest','off');        
        set(h_currentbar,'FaceColor', currentcolor, 'EdgeColor', 'none','HitTest','off');
        
        tidyAxes(metricWindow)
        set(metricWindow,'ytick',find(ismember(sorted_comps,current_comp)),'yticklabel','>')
        set(metricWindow,'xtick',[])
        set(metricWindow,'fontweight','bold','fontsize',16)
        axis(metricWindow,'tight')
        set(metricWindow,'ydir','reverse')      
        set(metricWindow,'ButtonDownFcn',@mouse_select_component)

        drawnow

        % Add info about metrics as title above tICWindow
        titlestr = 'Component ranking: ';
        for j = 1:length(metric_names)
            v = metrics.(metric_names{j}).value(current_comp); % Metric value
            r = sum(v <= metrics.(metric_names{j}).value);
            titlestr = [titlestr sprintf('%s: %i (%.2f)  ',metric_names{j},r,v)];
        end  
        title(tICWindow,titlestr,'fontsize',FONTSIZE,'interpreter','none')
        
    end


    function mouse_select_component(a,b,c)
        C = get(a,'CurrentPoint');
        current_comp_ind = round(C(1,2));
        current_comp = sorted_comps(current_comp_ind);
        redraw
    end


    
    function key_press(~,evnt)
        % scrolling through components
        if any(strcmp(evnt.Key,{'rightarrow','leftarrow','downarrow','uparrow'}))
            
            if strcmp(evnt.Key,'rightarrow') || strcmp(evnt.Key,'downarrow')
                current_comp_ind = current_comp_ind + 1;
                
            elseif strcmp(evnt.Key,'leftarrow') || strcmp(evnt.Key,'uparrow')
                current_comp_ind = current_comp_ind - 1;
                
            end
            
            if current_comp_ind < 1
                current_comp_ind = length(sorted_comps);
            elseif current_comp_ind > length(sorted_comps)
                current_comp_ind = 1;
            end

            current_comp = sorted_comps(current_comp_ind);
            redraw
        
        % hotkeys
        elseif strcmp(evnt.Key, 'a')
            current_metric = current_metric - 1;
            if current_metric <= 0 
                current_metric = length(metric_names);
            end
            set(uitools.metrics,'Value',current_metric);
            switchmetric
            
        elseif strcmp(evnt.Key, 's')
            current_metric = current_metric + 1;
            if current_metric > length(metric_names)
                current_metric = 1;
            end 
            set(uitools.metrics,'Value',current_metric);
            switchmetric
            
        elseif strcmp(evnt.Key, 'b')
            setbad
            
        elseif strcmp(evnt.Key, 'g')
            setgood
            
        elseif strcmp(evnt.Key, 'z')
            %cb_zoom - AB disabled until it can be made to set the button
            %toggle status too
            
        else
            % no action on other keys            
        end%if
    end % key_press



    function setgood(~,~)
        bad_components(bad_components==current_comp) = [];
        redraw
    end 

    function setbad(~,~)
        bad_components = unique(sort([bad_components(:); current_comp]));
        redraw
    end 


    function switchmetric(~,~)
        current_metric = get(uitools.metrics,'Value');
        [~,sorted_comps] = sort(metrics.(metric_names{current_metric}).value,'descend');
        
        switch lower(get(metricContext.switch,'checked'))
            case 'on'
                current_comp_ind = 1;
            case 'off'
                current_comp_ind = find(sorted_comps == current_comp); 
        end
        current_comp = sorted_comps(current_comp_ind);
        redraw
    end


    function cb_metricContext(src,~)
        
        switch lower(get(metricContext.switch,'checked'))
            case 'on'
                set(metricContext.switch,'checked','off');
            case 'off'
                set(metricContext.switch,'checked','on');
        end
    end


    function cb_zoom(~,~)
        zoom
    end


end


%% Really lazy hard coding of icons:
function [icon_good,icon_bad,icon_zoom] = load_icons


cdata = [
     0     0     0     0     0     0     0     0     0     0     0     1     1     1     1     1
     0     0     0     0     0     0     0     0     0     0     1     1     1     1     1     0
     0     0     0     0     0     0     0     0     0     1     1     1     1     1     0     0
     0     0     0     0     0     0     0     0     0     1     1     1     1     0     0     0
     0     0     0     0     0     0     0     0     1     1     1     1     1     0     0     0
     0     0     0     0     0     0     0     0     1     1     1     1     0     0     0     0
     0     0     0     0     0     0     0     1     1     1     1     1     0     0     0     0
     0     1     1     1     1     0     0     1     1     1     1     0     0     0     0     0
     0     1     1     1     1     1     1     1     1     1     1     0     0     0     0     0
     0     1     1     1     1     1     1     1     1     1     0     0     0     0     0     0
     0     0     1     1     1     1     1     1     1     1     0     0     0     0     0     0
     0     0     0     1     1     1     1     1     1     1     0     0     0     0     0     0
     0     0     0     0     1     1     1     1     1     0     0     0     0     0     0     0
     0     0     0     0     1     1     1     1     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0];
 
cdata_r = double(~cdata); cdata_r(cdata_r==1) = nan;
cdata_g = 0.5*double(~isnan(cdata));
cdata_b = double(~cdata);
icon_good = cat(3,cdata_r,cdata_g,cdata_b);

cdata = [
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     1     1     0     0     0     0     0     1     1     1     1     0
     0     0     1     1     1     1     1     0     0     0     1     1     1     1     1     0
     0     0     1     1     1     1     1     1     1     1     1     1     1     0     0     0
     0     0     0     1     1     1     1     1     1     1     1     1     0     0     0     0
     0     0     0     0     0     1     1     1     1     1     0     0     0     0     0     0
     0     0     0     0     1     1     1     1     1     1     1     0     0     0     0     0
     0     0     0     0     1     1     1     1     1     1     1     1     0     0     0     0
     0     0     0     1     1     1     1     1     0     0     1     1     1     0     0     0
     0     0     1     1     1     1     1     1     0     0     0     1     1     1     0     0
     0     0     1     1     1     1     1     0     0     0     0     0     1     1     1     0
     0     0     1     1     1     1     0     0     0     0     0     0     0     1     1     0
     0     0     1     1     1     1     0     0     0     0     0     0     0     1     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];

cdata_r = 0.8*cdata; cdata_r(cdata_r==0) = nan;
cdata_g = double(~cdata);
cdata_b = double(~cdata);
icon_bad = cat(3,cdata_r,cdata_g,cdata_b);
 
icon_zoom = load(fullfile(matlabroot,'toolbox','matlab','icons','zoom.mat'));
icon_zoom = icon_zoom.(char(fieldnames(icon_zoom)));

end

function [] = tidyAxes(h, FONTSIZE)
if nargin < 2, FONTSIZE = 14; end

set(h,...
    'FontName', 'Helvetica', ...
    'FontSize', FONTSIZE, ...
    'Box', 'on', ...
    'YGrid', 'on', ...
    'XGrid', 'off', ...
    'TickDir', 'in', ...
    'TickLength', [0.005 0.005], ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'XColor', [0.3 0.3 0.3], ...
    'YColor', [0.3 0.3 0.3], ...
    'LineWidth', 2);
end
