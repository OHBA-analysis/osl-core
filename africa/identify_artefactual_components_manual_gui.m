function bad_components = identify_artefactual_components_manual_gui(D,tc,topos,metrics,bad_components)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create UI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MainFig = figure('Name',            ['AFRICA: ' D.fname], ...
                 'NumberTitle',     'off'               , ...
                 'Menubar',         'none'              , ...
                 'DockControls',    'on'                , ...
                 'KeyPressFcn',     @key_press          , ...
                 'BusyAction',      'cancel'            , ...
                 'Visible',         'off'               , ...
                 'ResizeFcn',       @createlayout); % ,...
%'CloseRequestFcn', @close_fig);

goodcolor    = [ 48   128    20] / 255; % [0 0.5 0];50    205    50
badcolor     = [  0.8   0     0];
currentcolor = [100   149   237] / 255;  % [  0     0.9   1];
% tcplotcolor  = [ 39	   64   139] / 255;
pwplotcolor  = [100   149   237] / 255;
covcolor     = [121   121   121] / 255;
FONTSIZE     = 14;

%bad_components = ismember(1:size(tc,1),bad_components);

metric_names = fieldnames(metrics);
current_metric = 1;

[~,sorted_comps] = sort(metrics.(metric_names{current_metric}).value,'descend');

current_comp_ind = 1;
current_comp     = sorted_comps(current_comp_ind);


% Create plotting windows
tICWindow    = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');
specWindow   = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');
metricWindow = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');
covWindow    = axes('parent',MainFig, 'units', 'pixels','DrawMode','fast');
topoWindow = [];
for m = 1:size(topos,2)
    topoWindow(m) = axes('parent',MainFig, 'units','pixels', 'DrawMode','fast');
end

% icons
[icon_good,icon_bad,icon_zoom] = load_icons;

uitools.toolbar    = uitoolbar;
uitools.setgood    = uipushtool(uitools.toolbar,    'ClickedCallback',@setgood, 'CData',icon_good, 'TooltipString','Set component as good');
uitools.setbad     = uipushtool(uitools.toolbar,    'ClickedCallback',@setbad,  'CData',icon_bad,  'TooltipString','Set component as bad');    
uitools.zoom       = uitoggletool(uitools.toolbar,  'ClickedCallback',@cb_zoom, 'CData',icon_zoom,  'TooltipString','Set component as bad');    
uitools.metrics = uicontrol('Style', 'popup', 'String', metric_names, 'Position', [1 1 1 1], 'Callback', @switchmetric);   



% Add dropdown list of metrics to the toolbar
% warning off
% jToolbar = get(get(uitools.toolbar,'JavaContainer'),'ComponentPeer');
% if ~isempty(jToolbar)
%    jCombo = javax.swing.JComboBox(metric_names);
%    set(jCombo, 'ActionPerformedCallback', @switchmetric);
%    jToolbar(1).add(jCombo); 
%    jToolbar(1).repaint;
%    jToolbar(1).revalidate;
% end
% warning on



drawnow

redraw

jFig = get(handle(MainFig),'JavaFrame');
jFig.setMaximized(true);

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
        
    end




    function redraw % UPDATE PLOTS
        
        
        t = (1:size(tc,2))./D.fsample;
        
        if ismember(current_comp,bad_components)
            thistcplotcolor = badcolor;
            thispwplotcolor = badcolor;
        else
            thistcplotcolor = goodcolor; %tcplotcolor;
            thispwplotcolor = pwplotcolor;
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
            
            %colorbar('peer', topoWindow(mWindow)); % AB disabled due to
            %overlapping bars & topos for Elekta data

            
            %freezeColors

        end
        
        % Redraw spectrum window
        axes(specWindow); 
        pwelch(tc(current_comp,~isnan(tc(current_comp,:))),[],[],[],D.fsample); 
        set(findobj(specWindow,'type','line'),'color', thispwplotcolor);
        tidyAxes(specWindow, FONTSIZE);
        title('Component Power Spectrum', 'FontSize', FONTSIZE);
        
        % redraw covariate window
        axes(covWindow);
        isValidCovField =   isfield(metrics.(metric_names{current_metric}), 'timeCourse') && ...
                           ~isempty(metrics.(metric_names{current_metric}).timeCourse)    && ...
                             length(metrics.(metric_names{current_metric}).timeCourse) == length(t);
        if isValidCovField
            plot(t(:), metrics.(metric_names{current_metric}).timeCourse,'color',covcolor);
        else
            plot(t, zeros(size(t)), 'color', covcolor);
        end%if
        axis tight
        set(covWindow,'xlim',[t(1) t(end)]);
        set(covWindow,'ylim',max(abs(get(covWindow,'ylim')))*[-1 1]);
        title('Covariate Time Course', 'FontSize', FONTSIZE);
        xlabel('Time (s)', 'FontSize', FONTSIZE);
        tidyAxes(covWindow, FONTSIZE)
        
        % Redraw tIC window
        axes(tICWindow);  
        plot(t,tc(current_comp,:),'color',thistcplotcolor), 
        set(tICWindow,'ylim',max(abs(get(tICWindow,'ylim')))*[-1 1]);
        set(tICWindow,'xlim',[t(1) t(end)]);
        tidyAxes(tICWindow, FONTSIZE);
        title('Component Time Course', 'FontSize', FONTSIZE);
        %xlabel('Time (s)', 'FontSize', FONTSIZE);

        
        % Redraw metric window
        axes(metricWindow), cla(metricWindow), hold(metricWindow,'on');
        goodbars = metrics.(metric_names{current_metric}).value(sorted_comps); 
        goodbars(ismember(sorted_comps,bad_components))=nan;
        goodbars(sorted_comps==current_comp) = nan; 
        goodbars = barh(goodbars); set(goodbars, 'FaceColor', goodcolor) 
        badbars  = metrics.(metric_names{current_metric}).value(sorted_comps); 
        badbars(~ismember(sorted_comps,bad_components))=nan;
        badbars(sorted_comps==current_comp)=nan; 
        badbars  = barh(badbars);  set(badbars,  'FaceColor', badcolor) 
        currentbar = nan(size(sorted_comps));
        currentbar(sorted_comps==current_comp) = metrics.(metric_names{current_metric}).value(current_comp); 
        currentbar = barh(currentbar);  set(currentbar, 'FaceColor', currentcolor) 
        axis(metricWindow,'tight','off')
        set(metricWindow,'ydir','reverse')
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
        %current_metric = find(strcmp(metric_names,get(hCombo,'SelectedItem')));
        current_metric = get(uitools.metrics,'Value');
        [~,sorted_comps] = sort(metrics.(metric_names{current_metric}).value,'descend');
        current_comp = sorted_comps(1);
        current_comp_ind = 1;
        redraw
    end


    function cb_zoom(~,~)
        zoom
    end


end


%% Really lazy hard coding of icons:
function [icon_good,icon_bad,icon_zoom] = load_icons
clc

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
 
icon_zoom = load(fullfile(matlabroot,'/toolbox/matlab/icons/zoom.mat'));
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
