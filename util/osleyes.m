classdef osleyes < handle
    % A pure Matlab NIFTI viewer based on common OSL usage of fsleyes/fslview
    %
    % 
    % TERMINOLOGY
    % - A *layer* corresponds to a NIFTI file
    % - A *volume* corresponds to a slice in the 4th dimension of the data in the NIFTI file
    %
    % osleyes.colormaps - A cell array of colormaps for each layer. Each
    %                     element corresponds to a layer. The layer can
    %                     contain a colormap, or a cell array with 2 colormaps
    %                     for a bidirectional map. A colormap is 
    %                       - the name of a function that takes in a number of elements e.g. 'hot'
    %                       - a (N x 3) matrix of RGB values (e.g. the matrix returned by 'osl_colormap('hot')'
    %
    % OTHER NOTES
    % - figure contains 'osleyes' property with handle to this object
    % - figure will be closed if handle is deleted
    % - object will be deleted if figure is closed

    % Romesh Abeysuriya 2017

    % Properties controlling the GUI that users can interact with
    properties
        layer = struct() % Struct array with layer properties - colormap, clim, alpha, name, visible, volume
        current_point = [1 1 1]; % These are the XYZ MNI coordinates for the crosshair
        active_layer = 1; % The active layer controls which layer's properties are shown in the control panel and the colorbars
        show_controls = 1;
        show_crosshair = 1;
        title = ''; % String of title for plot
    end

    % Dependent properties the user might want to use in their code
    properties(Dependent)
        nvols % Number of volumes in the active layer
    end

    % Read-only properties the user might want to use in their code
    properties(SetAccess=protected)
        fig % Handle to the main figure window bound to the object
        images = {} % Input image data - file name, matrix, or struct
    end

    properties(GetAccess=private,SetAccess=private)
        img = struct('vol',{},'res',{},'xform',{},'toffset',{},'tunits',{}); % The image data

        ax % Handles of the three display axes
        controls % Handles for control panel and associated controls
        contextmenu % Handle for the context menu
        lims = nan(3,2); % Axis limits for each MNI dimension (used in plots)

        colormap_resolution = 255; % Number of colors in each colormap (if not specified as matrix)
        coord = {} % Axis coordinates for each image

        h_img = {} % Cell array of img handles associated with each nii file (there are 3)
        h_crosshair % Handles for crosshairs on each axis
        h_coloraxes % Handles to colorbar axes
        h_colorimage % Handles to colorbar images
        h_title % Handle to title uicontrol text object

        ts_ax  % Handle to timeseries axis
        ts_line  % Handle to line in timeseries plot
        ts_bar  % Handle to marker bar in timeseries plot
        ts_warning % Text object saying only one volume present
        ts_xlabel % Handle to xlabel on time series plot
        ts_title % Handle to timeseries title

        colormap_matrices = {}; % Cached colormaps - could be at higher resolution than input colormap
        under_construction = true; % Don't render anything while this is true
        motion_active = 0; % If this flag is true, then the slices will be updated when the mouse is moved

    end

    properties(GetAccess=private,SetAccess=private,Dependent)
        time % time values for the active layer
    end

    methods

        function self = osleyes(images,varargin)
            % Constructor for osleyes
            % 
            % INPUTS
            %
            % - images : specify input images. These can be in one of 4 formats
            %
            %   - A file name e.g. "myfile.nii.gz". A standard brain will
            %     automatically be guessed
            %   - A cell array of file names e.g.
            %     {"layer1.nii.gz","layer2.nii.gz"}. If the first entry is
            %     empty, then a standard mask/brain image will automatically
            %     be guessed.
            %   - A matrix can be used instead of strings to specify
            %     an image. In this case, the size of the matrix will be used
            %     to guess the standard mask. If no mask is found, an error
            %     will be raised. The matrix will automatically be reshaped if possible
            %   - A struct can be used instead of strings, containing fields
            %     'img' and 'res', and 'xform'. Can optionally also contain 'label'
            %
            %   This aims to provide a high degree of flexibility. For example, you could use
            %   
            %   osleyes({[],'myfile.nii.gz',x,struct('img',y,'xform',y_xform)})
            %   
            %   to mix all possible input types. Note that if you want to display an image
            % *without* a standard brain being automatically included, you need to input a
            % cell array e.g.
            %
            %   osleyes({[],'myfile.nii.gz'}) - display with standard brain as first layer
            %   osleyes('myfile.nii.gz') - exactly equivalent to above command
            %   osleyes({'myfile.nii.gz'}) - display only included file, no standard brain
            %
            % - varargin : This is a set of key-value pair
            %   pairs that are assigned to the layer properties or object after
            %   construction. For example osleyes('myfile','clims',{[],[1 1]})
            %
            % OUTPUTS
            % - instance of osleyes class

            if nargin < 1 || isempty(images) 
                images = {fullfile(osldir,'std_masks/MNI152_T1_2mm_brain.nii.gz')};
            end
            
            if ~iscell(images)
                images = {[],images};
            end

            if isempty(images{1})
                tmp_img = load_image(images{2});
                try
                    [~,images{1},~] = parcellation.guess_template(tmp_img.vol);
                catch ME
                    switch ME.identifier
                        case 'osl:parcellation:no_matching_mask'
                            fprintf(2,'Warning - No matching standard brain mask was found\n');
                            images = images(2:end);
                        otherwise
                            rethrow(ME)
                    end
                end
            end

            self.fig = figure('Units','Characters','Color','k','Menubar','none','InvertHardCopy','off');
            self.initial_render();
            set(self.fig,'KeyPressFcn',@(a,b) KeyPressFcn(self,a,b),'CloseRequestFcn',@(~,~) delete(self),'ResizeFcn',@(~,~) resize(self));
            addprop(self.fig,'osleyes');
            set(self.fig,'osleyes',self); % Store handle to this osleyes in the figure so it can be retrieved later if desired

            self.images = images;
            self.controls.image_list.String = {};

            for j = 1:length(self.images)

                [self.img(j),self.layer(j).name] = load_image(self.images{j});
                self.coord{j} = get_coords(self.img(j).vol,self.img(j).xform);
                self.h_img{j}(1) = image(self.coord{j}.y,self.coord{j}.z,permute(self.img(j).vol(1,:,:,1),[2 3 1])','Parent',self.ax(1),'HitTest','off','AlphaDataMapping','none','AlphaData',1);
                self.h_img{j}(2) = image(self.coord{j}.x,self.coord{j}.z,permute(self.img(j).vol(:,1,:,1),[1 3 2])','Parent',self.ax(2),'HitTest','off','AlphaDataMapping','none','AlphaData',1);
                self.h_img{j}(3) = image(self.coord{j}.x,self.coord{j}.y,permute(self.img(j).vol(:,:,1,1),[1 2 3])','Parent',self.ax(3),'HitTest','off','AlphaDataMapping','none','AlphaData',1);
                
                self.layer(j).clim = [0,max(self.img(j).vol(:))]; % Default color range

                if min(self.img(j).vol(:)) < 0
                    self.layer(j).colormap = {osl_colormap('hot'),osl_colormap('cold')};
                elseif j == 1
                    self.layer(j).colormap = osl_colormap('grey');
                else
                    self.layer(j).colormap = osl_colormap('hot');
                end

                if isempty(self.layer(j).name)
                    self.layer(j).name = sprintf('Image %d',j);
                end

                self.layer(j).volume = 1;
                self.layer(j).visible = 1;
                self.layer(j).alpha = 1;
            end

            self.lims(1,:) = [min(cellfun(@(x) min(x.x),self.coord)) max(cellfun(@(x) max(x.x),self.coord))];
            self.lims(2,:) = [min(cellfun(@(x) min(x.y),self.coord)) max(cellfun(@(x) max(x.y),self.coord))];
            self.lims(3,:) = [min(cellfun(@(x) min(x.z),self.coord)) max(cellfun(@(x) max(x.z),self.coord))];
            set(self.ax(1),'XLim',self.lims(2,:),'YLim',self.lims(3,:),'Visible','off','Color','k','Clipping','off','View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal');
            set(self.ax(2),'XLim',self.lims(1,:),'YLim',self.lims(3,:),'Visible','off','Color','k','Clipping','off','View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
            set(self.ax(3),'XLim',self.lims(1,:),'YLim',self.lims(2,:),'Visible','off','Color','k','Clipping','off','View', [0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
            orientation_letters(self.ax(1),{'P','A','I','S'});
            orientation_letters(self.ax(2),{'L','R','I','S'});
            orientation_letters(self.ax(3),{'L','R','P','A'});

            self.h_crosshair(1) = plot(self.ax(1),NaN,NaN,'g','HitTest','off');
            self.h_crosshair(2) = plot(self.ax(2),NaN,NaN,'g','HitTest','off');
            self.h_crosshair(3) = plot(self.ax(3),NaN,NaN,'g','HitTest','off');

            set(self.fig,'WindowButtonDownFcn',@(~,~) self.activate_motion())
            set(self.fig,'WindowButtonMotionFcn',@(a,b) move_mouse(self))
            set(self.fig,'WindowButtonUpFcn',@(~,~) self.deactivate_motion())

            self.contextmenu.root = uicontextmenu();
            self.contextmenu.plot_timeseries = uimenu(self.contextmenu.root, 'label','Plot time series','Enable','On','Callback',@(~,~) self.plot_timeseries());
            set(self.fig,'uicontextmenu',self.contextmenu.root);

            self.resize()
            self.under_construction = false; % Enable rendering
            self.active_layer = length(self.images); % Triggers re-rendering

            % Set options
            arg = inputParser;
            arg.KeepUnmatched = true;
            arg.parse(varargin{:});
            f = fields(arg.Unmatched);
            for j = 1:length(f)
                if ismember(f{j},fieldnames(self.layer))
                    for k = 1:length(self.layer)
                        if iscell(arg.Unmatched.(f{j}))
                            if ~isempty(arg.Unmatched.(f{j}){k})
                                self.layer(k).(f{j}) = arg.Unmatched.(f{j}){k};
                            end
                        else
                            if isfinite(arg.Unmatched.(f{j})(k)) 
                                self.layer(k).(f{j}) = arg.Unmatched.(f{j})(k);
                            end
                        end
                    end
                else
                    self.(f{j}) = arg.Unmatched.(f{j});
                end
            end
        end

        function delete(self)
            % Destructor closes the figure when the object ceases to exist
            delete(self.fig);
        end

        function animate(self,fps)
            % Cycle through volumes
            % Pause assuming no rendering time - so fps is only approximate
            if nargin < 2 || isempty(fps) 
                fps = 30;
            end
            
            if self.nvols == 1
                error('Layer only has one volume, there is nothing to animate!');
            end

            while 1
                self.layer(self.active_layer).volume = 1+mod(self.layer(self.active_layer).volume,self.nvols);
                pause(1/fps);
            end
        end
        
        function plot_timeseries(self)
            % This function will open a new figure window
            % Plotting the timeseries itself is handled in two parts
            % - The line is refreshed by refresh_slices when the MNI coordinates change
            % - The bar is refreshed by set.current_vols when the volume changes

            % Odd things will happen if more than one window is allowed to be bound to the
            % same osleyes instance
            if ishandle(self.ts_line)
                return
            end

            fig = figure;
            self.ts_ax = axes(fig);
            self.ts_line = plot(self.ts_ax,1,1,'HitTest','off');
            hold(self.ts_ax,'on')
            self.ts_bar = plot(self.ts_ax,1,1,'r','HitTest','off');
            self.ts_warning = uicontrol(fig,'style','text','String','ONLY ONE VOLUME PRESENT IN ACTIVE LAYER','Units','normalized','Position',[0.25 0.25 0.5 0.5],'HitTest','off','FontSize',20,'BackgroundColor','w');
            set(self.ts_ax,'ButtonDownFcn',@(~,~) self.ts_set_volume());
            self.ts_xlabel = xlabel(self.ts_ax,'Volume/Time');
            self.ts_title = title(self.ts_ax,'Timeseries','Interpreter','none');
            self.active_layer = self.active_layer; % Reset the data in h_bar via set.current_vols and set axis limits

        end

        function set.layer(self,val)
            % Generic layer change function when clim, colormap, visible, alpha, vol are changed
            % UPDATE HIERARCHY
            % - active_layer (includes colormaps, so is slower) 
            %   - refresh_slices (only re-renders the image, must be fast because called while dragging)

            if self.under_construction
                self.layer = val;
                return
            end

            assert(length(val)==length(self.img),'Number of layer entries does not match number of images')
            assert(isempty(setdiff(fieldnames(val),fieldnames(self.layer))),'New layer fields do not match - fields must be: %s', cell2mat(cellfun(@(x) [x,', '],fieldnames(self.layer),'UniformOutput',false)'));

            colors_required = false; % Flag to check if we need to recompute colormaps
            for j = 1:length(val)

                % Check colormaps
                assert(~iscell(val(j).colormap) || (iscell(val(j).colormap) && length(val(j).colormap)==2),'Colormap must either be a value or a cell array of length 2');
                if ~strcmp(class(val(j).colormap),class(self.layer(j).colormap))
                    colors_required = true;
                elseif iscell(val(j).colormap)
                    colors_required = colors_required || any(size(val(j).colormap{1})~=size(self.layer(j).colormap{1})) || any(size(val(j).colormap{2})~=size(self.layer(j).colormap{2})) || any(val(j).colormap{1}(:)~=self.layer(j).colormap{1}(:)) || any(val(j).colormap{2}(:)~=self.layer(j).colormap{2}(:));
                else
                    colors_required = colors_required || any(size(val(j).colormap)~=size(self.layer(j).colormap)) || any(val(j).colormap(:)~=self.layer(j).colormap(:));
                end

                % Check clims
                assert(length(val(j).clim)==2,'Colour limits must have two elements')
                assert(val(j).clim(2)>=val(j).clim(1),'Colour limits must be in ascending order')
                colors_required = colors_required || any(val(j).clim(:)~=self.layer(j).clim(:));

                % Check alpha
                assert(isscalar(val(j).alpha) && isnumeric(val(j).alpha) && isfinite(val(j).alpha))
                assert(val(j).alpha>=0 && val(j).alpha<=1)

                % Check volume indices
                assert(isscalar(val(j).volume) && isnumeric(val(j).volume) && isfinite(val(j).volume))
                assert(val(j).volume>0 && val(j).volume <= size(self.img(j).vol,4),'Volume index out of bounds - this layer only has %d volumes',size(self.img(j).vol,4));

                % Check name
                assert(ischar(val(j).name) || isstring(val(j).name))
                colors_required = colors_required || ~strcmp(val(j).name,self.layer(j).name);

                % Check visible
                assert(isscalar(val(j).visible))
            end

            self.layer = val;

            if colors_required
                self.active_layer = self.active_layer;
            else
                self.refresh_slices();
            end

        end

        function set.current_point(self,val)
            % Validate limits
            xl = get(self.ax(3),'XLim');
            yl = get(self.ax(1),'XLim');
            zl = get(self.ax(1),'YLim');

            % Check that coordinates are within axis limits in all axes
            if val(1) > xl(2) || val(1) < xl(1) || val(2) > yl(2) || val(2) < yl(1) || val(3) > zl(2) || val(3) < zl(1)
                return
            end

            self.current_point = val;
            refresh_slices(self)
        end

        function v = set.active_layer(self,val)
            % Refresh the UI - called also via set.layer()
            assert(val > 0,'Layer must be positive');
            assert(val <= length(self.images),'Layer number cannot exceed number of layers');
            assert(val == round(val),'Layer must be an integer');

            % Refresh the colormaps
            for j = 1:length(self.layer)
                self.colormap_matrices{j} = self.compute_color_matrix(self.layer(j).colormap);
            end

            self.active_layer = val;
            set(self.controls.image_list,'Value',val,'String',{self.layer.name})
            set(self.controls.clim(1),'String',sprintf('%g',self.layer(self.active_layer).clim(1)));
            set(self.controls.clim(2),'String',sprintf('%g',self.layer(self.active_layer).clim(2)));
            set(self.controls.volume,'String',sprintf('%d',self.layer(self.active_layer).volume));
            set(self.controls.volume_label_count,'String',sprintf('of %d',self.nvols));

            set(self.controls.visible,'Value',self.layer(self.active_layer).visible);

            tickstrs = @(low,high,n) arrayfun(@(x) sprintf('%.2g',x),linspace(low,high,n),'UniformOutput',false);

            if iscell(self.layer(self.active_layer).colormap)
                set(self.h_coloraxes(1),'Visible','on');
                set(self.h_colorimage(1),'Visible','on');
                set(self.h_colorimage(2),'CData',permute(self.colormap_matrices{self.active_layer}{1},[1 3 2]));
                set(self.h_colorimage(1),'CData',permute(self.colormap_matrices{self.active_layer}{2},[1 3 2]));
                set(self.h_coloraxes(1),'YTick',linspace(0,1,4),'YTickLabel',tickstrs(-self.layer(self.active_layer).clim(1),-self.layer(self.active_layer).clim(2),4))
                set(self.h_coloraxes(2),'YTick',linspace(0,1,4),'YTickLabel',tickstrs(self.layer(self.active_layer).clim(1),self.layer(self.active_layer).clim(2),4))
            else
                set(self.h_coloraxes(1),'Visible','off');
                set(self.h_colorimage(1),'Visible','off');
                set(self.h_colorimage(2),'CData',permute(self.colormap_matrices{self.active_layer},[1 3 2]));
                set(self.h_coloraxes(2),'YTick',linspace(0,1,4),'YTickLabel',tickstrs(self.layer(self.active_layer).clim(1),self.layer(self.active_layer).clim(2),4))
            end

            if self.nvols==1
                set(self.controls.volume,'Enable','off');
            else
                set(self.controls.volume,'Enable','on');
            end

            % Hide time if time units are not defined
            if isempty(self.img(self.active_layer).tunits)
                set(self.controls.marker(5),'Visible','off');
            else
                set(self.controls.marker(5),'Visible','on');
            end

            if ishandle(self.ts_line)
                set(self.ts_line,'XData',self.time,'YData',nan(1,self.nvols));
                set(self.ts_ax,'XLim',[self.time(1) self.time(end)+1e-5*(self.time(1)==self.time(end))],'YLim',[min(self.img(self.active_layer).vol(:)) max(self.img(self.active_layer).vol(:)) ]);
            end

            self.resize(); % Update sizes of colorbars
            self.refresh_slices();
        end

        function set.title(self,str)
            set(self.h_title,'String',str);
            self.resize() % Resize the GUI in case the visibility of the title bar has changed
        end

        function set.show_controls(self,val)
            self.show_controls = logical(val);
            self.resize()
        end

        function set.show_crosshair(self,val)
            self.show_crosshair = logical(val);
            if self.show_crosshair
                set(self.h_crosshair,'Visible','on');
            else
                set(self.h_crosshair,'Visible','off');
            end
            self.refresh_slices()
        end

        function n = get.nvols(self)
            n = size(self.img(self.active_layer).vol,4);
        end

        function t = get.time(self)
            t = (0:self.nvols-1)*self.img(self.active_layer).res(4)+self.img(self.active_layer).toffset;
        end

    end

    methods(Access=private)

        function refresh_slices(self)
            % Update the slices without re-rendering colour bars
            % This will move the crosshairs, set the position markers, and change which slice is displayed

            p = self.current_point; % Current point in 3D

            set(self.h_crosshair(1),'XData',[self.lims(2,:) NaN p(2) p(2)],'YData',[p(3) p(3) NaN self.lims(3,:)]);
            set(self.h_crosshair(2),'XData',[self.lims(1,:) NaN p(1) p(1)],'YData',[p(3) p(3) NaN self.lims(3,:)]);
            set(self.h_crosshair(3),'XData',[self.lims(1,:) NaN p(1) p(1)],'YData',[p(2) p(2) NaN self.lims(2,:)]);
            set(self.controls.marker(1),'String',sprintf('X = %+06.1f',p(1)));
            set(self.controls.marker(2),'String',sprintf('Y = %+06.1f',p(2)));
            set(self.controls.marker(3),'String',sprintf('Z = %+06.1f',p(3)));
            set(self.controls.marker(5),'String',sprintf('Time = %+03.3f%s',self.img(self.active_layer).toffset+(self.layer(self.active_layer).volume-1)*self.img(self.active_layer).res(4),self.img(self.active_layer).tunits)); 

            set(self.controls.volume,'String',sprintf('%d',self.layer(self.active_layer).volume));

            % Now update each slice
            for j = 1:length(self.img)
                [~,idx(1)] = min(abs(self.coord{j}.x-p(1)));
                [~,idx(2)] = min(abs(self.coord{j}.y-p(2)));
                [~,idx(3)] = min(abs(self.coord{j}.z-p(3)));

                % These are the slice maps - need to convert them to color values now
                d1 = permute(self.img(j).vol(idx(1),:,:,self.layer(j).volume),[2 3 1])';
                d2 = permute(self.img(j).vol(:,idx(2),:,self.layer(j).volume),[1 3 2])';
                d3 = permute(self.img(j).vol(:,:,idx(3),self.layer(j).volume),[1 2 3])';

                if j == self.active_layer
                    set(self.controls.marker(4),'String',sprintf('Value = %+ 7.2f',self.img(j).vol(idx(1),idx(2),idx(3),self.layer(j).volume)));
                end

                % If bidirectional colormap, then hide abs(x)<clim(1), otherwise hide x<clim(1)
                if iscell(self.layer(j).colormap)
                    hidefcn = @abs; 
                else
                    hidefcn = @(x) x;
                end

                if self.layer(j).visible
                    set(self.h_img{j}(1),'Visible','on','CData',map_colors(d1,self.colormap_matrices{j},self.layer(j).clim),'AlphaData',self.layer(j).alpha*+(hidefcn(d1)>self.layer(j).clim(1)));
                    set(self.h_img{j}(2),'Visible','on','CData',map_colors(d2,self.colormap_matrices{j},self.layer(j).clim),'AlphaData',self.layer(j).alpha*+(hidefcn(d2)>self.layer(j).clim(1)));
                    set(self.h_img{j}(3),'Visible','on','CData',map_colors(d3,self.colormap_matrices{j},self.layer(j).clim),'AlphaData',self.layer(j).alpha*+(hidefcn(d3)>self.layer(j).clim(1)));
                else
                    set(self.h_img{j}(1),'Visible','off');
                    set(self.h_img{j}(2),'Visible','off');
                    set(self.h_img{j}(3),'Visible','off');
                end
            end

            if ishandle(self.ts_line)
                [~,idx(1)] = min(abs(self.coord{self.active_layer}.x-p(1)));
                [~,idx(2)] = min(abs(self.coord{self.active_layer}.y-p(2)));
                [~,idx(3)] = min(abs(self.coord{self.active_layer}.z-p(3)));
                set(self.ts_line,'YData',squeeze(self.img(self.active_layer).vol(idx(1),idx(2),idx(3),:)));
                set(self.ts_bar,'XData',[1 1]*self.time(self.layer(self.active_layer).volume),'YData',get(get(self.ts_bar,'Parent'),'YLim'));
                set(self.ts_title,'String',sprintf('%s - MNI (%.2f,%.2f,%.2f)',self.controls.image_list.String{self.active_layer},p(1),p(2),p(3)));

                if isempty(self.img(self.active_layer).tunits)
                    set(self.ts_xlabel,'String','Volume');
                else
                    set(self.ts_xlabel,'String',sprintf('Time (%s)',self.img(self.active_layer).tunits));
                end

                if self.nvols==1
                    set(self.ts_warning,'Visible','on');
                    set(self.ts_bar,'Visible','off')
                else
                    set(self.ts_warning,'Visible','off');
                    set(self.ts_bar,'Visible','on')
                end
            end
        end

        function ts_set_volume(self)
            % This callback runs when the user clicks on the timeseries to change the volume
            if isvalid(self) && self.nvols > 1
                p = get(self.ts_ax,'CurrentPoint');
                idx = (p(1,1)-self.img(self.active_layer).toffset)/self.img(self.active_layer).res(4);
                self.layer(self.active_layer).volume = ceil(idx);
            end

        end

        function activate_motion(self)  
            % This function is called whenever the user presses the mouse button
            % It sets the motion_active flag from 0 to the index of which axis the user
            % clicked on (and thus which axis should be used to compute the coordinates from)
            if strcmp(get(self.fig,'SelectionType'),'alt')
                return
            end

            if is_within(self.fig,self.ax(1))
                self.motion_active = 1;
            elseif is_within(self.fig,self.ax(2))
                self.motion_active = 2;
            elseif is_within(self.fig,self.ax(3))
                self.motion_active = 3;
            end
            move_mouse(self)
        end

        function deactivate_motion(self)
            % This function is called when the user releases the mouse button
            self.motion_active = 0;
        end

        function initial_render(self)
            % Create all of the objects AND perform the one-off layout of the control panel
            self.ax(1) = axes('Parent',self.fig,'Units','characters');
            self.ax(2) = axes('Parent',self.fig,'Units','characters');
            self.ax(3) = axes('Parent',self.fig,'Units','characters');
            self.controls.panel = uipanel(self.fig,'BorderType','none','Units','characters');
            hold(self.ax(1),'on');
            hold(self.ax(2),'on');
            hold(self.ax(3),'on');
            self.h_coloraxes(1) = axes('Box','on','Color','k','Units','characters');
            self.h_coloraxes(2) = axes('Box','on','Color','k','Units','characters');
            self.h_colorimage(1) = image([0 1],[0 1],1,'Parent',self.h_coloraxes(1));
            self.h_colorimage(2) = image([0 1],[0 1],1,'Parent',self.h_coloraxes(2));
            set(self.h_coloraxes,'XLim',[0 1],'YLim',[0 1],'XColor','w','YColor','w','XTick',[],'YDir','reverse','YAxisLocation','right');
            set(self.h_coloraxes(2),'YDir','normal');
            
            self.h_title = uicontrol(self.fig,'style','text','String','','FontSize',12,'Units','characters','Position',[0 0 0 0],'HorizontalAlignment','center','FontWeight','bold','BackgroundColor','k','ForegroundColor','w','HitTest','off');

            self.controls.image_list = uicontrol(self.controls.panel,'Callback',@(~,~) image_list_callback(self),'style','popupmenu','String','test','Units','characters','Position',[0 0.75 20 1.5]);

            self.controls.clim(1) = uicontrol(self.controls.panel,'Callback',@(~,~) clim_box_callback(self),'style','edit','String','1.0','Units','characters','Position',[0 0.2 7 1.2]);
            self.controls.clim(2) = uicontrol(self.controls.panel,'Callback',@(~,~) clim_box_callback(self),'style','edit','String','1.2','Units','characters','Position',[0 0.2 7 1.2]);
            self.controls.clim_label(1) = uicontrol(self.controls.panel,'style','text','String','Range:','Units','characters','Position',[0 0.3 8 1]);
            self.controls.clim_label(2) = uicontrol(self.controls.panel,'style','text','String','to','Units','characters','Position',[0 0.3 3 1],'HorizontalAlignment','left');

            self.controls.volume_label = uicontrol(self.controls.panel,'style','text','String','Volume:','Units','characters','Position',[0 1.6 8 1]);
            self.controls.volume = uicontrol(self.controls.panel,'Callback',@(~,~) volume_box_callback(self),'style','edit','String','1','Units','characters','Position',[0 1.5 7 1.2]);
            self.controls.volume_label_count = uicontrol(self.controls.panel,'style','text','String','of XXXX','Units','characters','Position',[0 1.6 8 1],'HorizontalAlignment','left');

            self.controls.visible = uicontrol(self.controls.panel,'Callback',@(~,~) visible_box_callback(self),'style','checkbox','Units','characters','Position',[0 1 3 1]);

            self.controls.marker(1) = uicontrol(self.controls.panel,'style','text','String','X = +000.0','Units','characters','Position',[0 2 11 1],'HorizontalAlignment','left');
            self.controls.marker(2) = uicontrol(self.controls.panel,'style','text','String','Y = +000.0','Units','characters','Position',[0 1 11 1],'HorizontalAlignment','left');
            self.controls.marker(3) = uicontrol(self.controls.panel,'style','text','String','Z = +000.0','Units','characters','Position',[0 0 11 1],'HorizontalAlignment','left');
            self.controls.marker(4) = uicontrol(self.controls.panel,'style','text','String','Value = +0000.0','Units','characters','Position',[0 1 18 1],'HorizontalAlignment','left');
            self.controls.marker(5) = uicontrol(self.controls.panel,'style','text','String','(Time = +000.0s)','Units','characters','Position',[0 2 18 1],'HorizontalAlignment','left','Visible','off'); % Start invisible to avoid flicker

            % Put b to the right of a with given padding
            next_to = @(b,a,padding)    set(b,'Position',sum(get(a,'Position').*[1 0 1 0]).*[1 0 0 0]  + [padding 0 0 0] + [0 1 1 1].*get(b,'Position'));

            set(self.controls.visible,'Position',[1 1 3 1])
            next_to(self.controls.image_list,self.controls.visible,1);

            next_to(self.controls.volume_label,self.controls.image_list,1);
            next_to(self.controls.volume,self.controls.volume_label,1);
            next_to(self.controls.volume_label_count,self.controls.volume,1);

            next_to(self.controls.clim_label(1),self.controls.image_list,1);
            next_to(self.controls.clim(1),self.controls.clim_label(1),1);
            next_to(self.controls.clim_label(2),self.controls.clim(1),1);
            next_to(self.controls.clim(2),self.controls.clim_label(2),0);

            next_to(self.controls.marker(1),self.controls.clim(2),2);
            next_to(self.controls.marker(2),self.controls.clim(2),2);
            next_to(self.controls.marker(3),self.controls.clim(2),2);
            next_to(self.controls.marker(4),self.controls.marker(3),0.5);
            next_to(self.controls.marker(5),self.controls.marker(3),0.5);

        end


        function resize(self)
            % Resize the main objects in the figure window
            % The contents of the control panel will not be changed - their position is determined
            % by the initial layout function and is fixed from that point on
            orig_units = get(self.fig,'Units');
            set(self.fig,'Units','characters');
            figpos = get(self.fig,'Position');

            % Layout is - 0.3 reserved for each axis, plus 0.1 for the two colorbars
            cb = 0.1; % Colorbar width          
            w = 0.3*(1-cb)*figpos(3);
            p = ((1-cb)*figpos(3)-3*w)/6;
            control_height = 3*+self.show_controls; % Control panel height in characters
            title_height = 2*~isempty(get(self.h_title,'String'));
            ax_height = figpos(4)-control_height-title_height;

            if ax_height <= 0 % If user made the window extremely small
                return
            end

            set(self.h_title,'Position',[0 figpos(4)-title_height figpos(3) title_height-0.5*(title_height~=0)]);
            set(self.ax(1),'Position',[p control_height w ax_height]);
            set(self.ax(2),'Position',[p+w+p+p  control_height w ax_height]);
            set(self.ax(3),'Position',[p+w+p+p+w+p+p control_height w ax_height]);
            set(self.controls.panel,'Position',[0 0 figpos(3) control_height]);

            if iscell(self.layer(self.active_layer).colormap)
                set(self.h_coloraxes(1),'Position',[(1-cb+0.005)*figpos(3) control_height+ax_height*0.075 0.015*figpos(3) ax_height*0.4]);
                set(self.h_coloraxes(2),'Position',[(1-cb+0.005)*figpos(3) control_height+ax_height*0.525 0.015*figpos(3) ax_height*0.4]);
            else
                set(self.h_coloraxes(2),'Position',[(1-cb+0.005)*figpos(3) control_height+ax_height*0.1 0.015*figpos(3) ax_height*0.8]);
            end
            set(self.fig,'Units',orig_units);
        end


        function m = compute_color_matrix(self,cmap)
            if iscell(cmap)
                m = cellfun(@(x) self.compute_color_matrix(x),cmap,'UniformOutput',false);
                return
            end

            if isstr(cmap)
                m = feval(cmap,self.colormap_resolution);
            else
                assert(size(cmap,2)==3,'Colormap must have 3 columns (RGB matrix)')
                m = cmap;
            end
        end

    end

end

function image_list_callback(self)
    % This callback is run when the user selects a layer from the dropdown list
    self.active_layer = get(self.controls.image_list,'Value');
end

function clim_box_callback(self)
    % This callback is run whenever EITHER of the colour range textboxes is modified
    self.layer(self.active_layer).clim = [str2double(get(self.controls.clim(1),'String')) str2double(get(self.controls.clim(2),'String'))];
end

function visible_box_callback(self)
    % This callback is run if the user clicks the visibility checkbox
    self.layer(self.active_layer).visible = get(self.controls.visible,'Value');
end

function volume_box_callback(self)
    % If the user types a umber into the volume edit text box, this callback is run
    new_volume = str2double(get(self.controls.volume,'String'));
    if new_volume < 1 || new_volume > self.nvols
        uiwait(errordlg(sprintf('Volume out of bounds (image has %d volumes)',self.nvols)));
        set(self.controls.volume,'String',num2str(self.layer(self.active_layer).volume));
    else
        self.layer(self.active_layer).volume = new_volume;
    end
end

function rgb = map_colors(x,cmap,clim)
    % Map the image to colours using the colour matrices
    % The problem is that each axis contains multiple images each with 
    % a different colourmap. But an axis can only have one colormap. Therefore
    % we need to perform the colour mapping ourselves
    % 
    % INPUT
    %   x - image data
    %   cmap - color matrix OR cell with two color matrices
    %   clim - color limit for this plot
    % OUTPUT
    %   rgb - colour image data for rendering to screen

    if iscell(cmap) % We need to treat the positive and negative parts separately
        % These are the positive indices
        pos_idx = min(size(cmap{1},1),round((size(cmap{1},1)-1)*(x-clim(1))/(clim(2)-clim(1)))+1);
        pos_idx(pos_idx<1) = 1;
        neg_idx = min(size(cmap{2},1),round((size(cmap{2},1)-1)*(-x-clim(1))/(clim(2)-clim(1)))+1);
        neg_idx(neg_idx<1) = 1;

        r = zeros(size(x));
        g = zeros(size(x));
        b = zeros(size(x));
        r(x>=clim(1)) = cmap{1}(pos_idx(x>=clim(1)),1);
        g(x>=clim(1)) = cmap{1}(pos_idx(x>=clim(1)),2);
        b(x>=clim(1)) = cmap{1}(pos_idx(x>=clim(1)),3);
        r(x<clim(1)) = cmap{2}(neg_idx(x<clim(1)),1);
        g(x<clim(1)) = cmap{2}(neg_idx(x<clim(1)),2);
        b(x<clim(1)) = cmap{2}(neg_idx(x<clim(1)),3);
    else
        idx = min(size(cmap,1),round((size(cmap,1)-1)*(x-clim(1))/(clim(2)-clim(1)))+1);
        idx(idx<1) = 1;

        r = reshape(cmap(idx,1),size(x));
        g = reshape(cmap(idx,2),size(x));
        b = reshape(cmap(idx,3),size(x));
    end

    rgb = cat(3,r,g,b);

end

function c = get_coords(img,xform)
    % Take in an image and xform matrix
    % Return a struct with fields 'x','y','z' containing vectors
    % of MNI coordinates for each dimension. These are then used by image() when rendering
    % the layer to get the spatial position, orientation, and size correct
    %
    % See https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html#ref0
    %
    % Note - this is the point at which xform is used to quantify how the
    % left-right axes are reversed The arrays here are used as the X and Y
    % data values in calls to imagesc() which ensures that the image is
    % correctly oriented relative to its parent axis. Then, the *axis* can be
    % flipped left or right as desired to display in radiological orientation
    % or not The essential thing is that here, the order of the arrays encodes
    % whether they are in ascending or descending order in MNI space

    i = 0:size(img,1)-1;
    j = 0:size(img,2)-1;
    k = 0:size(img,3)-1;

    assert(isdiag(xform(1:3,1:3)),'xform matrix is not diagonal which means there are complex transformations in the NIFTI file. OSLEYES cannot handle these transformations, you must use FSLEYES')
    assert(all(xform(4,:)==[0 0 0 1]),'Unexpected last row of xform - was expecting [0 0 0 1]. OSLEYES has not been tested with this type of file. Use FSLEYES')

    c.x = xform(1,1) * i + xform(1,4);
    c.y = xform(2,2) * j + xform(2,4);
    c.z = xform(3,3) * k + xform(3,4);

end

function move_mouse(self)
    % This function is called whenever the mouse is moved and is responsible for
    % redrawing the slices via set.current_point()

    if isMultipleCall();  return;  end
    if ~self.motion_active; return; end

    if self.motion_active == 1 && is_within(self.fig,self.ax(1))
        p = get(self.ax(1),'CurrentPoint');
        self.current_point(2:3) = p(1,1:2);
    elseif self.motion_active == 2 && is_within(self.fig,self.ax(2))
        p = get(self.ax(2),'CurrentPoint');
        self.current_point([1,3]) = p(1,1:2);
    elseif self.motion_active == 3 && is_within(self.fig,self.ax(3))
        p = get(self.ax(3),'CurrentPoint');
        self.current_point(1:2) = p(1,1:2);
    end
end

function within = is_within(f,h)
    % This function takes in 
    % - f - figure handle
    % - h - handle of object within f
    % and returns true if the CurrentPoint of the figure is within the 
    % bounds of the object h (e.g. if the mouse is over an axis)

    C = get(f,'CurrentPoint');
    pos = get(h,'Position');

    if C(1) < pos(1) || C(1) > pos(1)+pos(3) || C(2) < pos(2) || C(2) > pos(2)+pos(4)
        within = false;
    else
        within = true;
    end
end

function flag=isMultipleCall()
    % Based on Yair Altman's function
    flag = false; 
    % Get the stack
    s = dbstack();
    if numel(s)<=2
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

function orientation_letters(ax,labels)
    % This function adds the orientation letter marker to the main display
    %
    % Labels will be positioned at:
    % {left, right, bottom, top}
    % Hold will be left on
    xl = get(ax,'XLim');
    yl = get(ax,'YLim');
    hold(ax,'on')
    %text(min(xl),mean(yl),labels{1},'Parent',ax,'Color','w','FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle','HitTest','off','FontSize',10)
    text(max(xl),mean(yl),labels{2},'Parent',ax,'Color','w','FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','middle','HitTest','off','FontSize',10)
    text(mean(xl),min(yl),labels{3},'Parent',ax,'Color','w','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom','HitTest','off','FontSize',10)
    %text(mean(xl),max(yl),labels{4},'Parent',ax,'Color','w','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','HitTest','off','FontSize',10)
end


function KeyPressFcn(self,~,KeyData)
    switch KeyData.Key
        case 'uparrow'
            self.layer(self.active_layer).volume = 1+mod(self.layer(self.active_layer).volume,self.nvols);
        case 'downarrow'
            self.layer(self.active_layer).volume = 1+mod(self.layer(self.active_layer).volume-2,self.nvols);
    end
end

function [img,name] = load_image(x)
    % Load an input image, return img struct
    % Supported input formats are
    % - string - load NII file
    % - struct - populate fields img, xform (mandatory), name='', res=[1 1 1 1], tunits = '', toffset=0 (optional),
    % - matrix - try and guess the standard mask and use that to retrieve xform etc.
    %
    % This function should produce sensible outputs except
    % if the input image is a 2D matrix whose row count does not match
    % any standard mask - in which case, turning into a 3D matrix 
    % is impossible and a hard error will be thrown
    name = '';

    img = struct();

    if ischar(x) || isstring(x) % Load a NIFTI file
        [img.vol,img.res,img.xform,~,img.toffset,img.tunits] = nii.load(x);
        [~,fname,ext] = fileparts(x);
        if isempty(ext)
            name = fname;
        else
            name = [fname '.' ext];
        end
    elseif isstruct(x)
        img.vol = x.img;
        img.xform = x.xform;

        if ~isfield(x,'res')
            img.res = [1 1 1 1];
        else
            img.res = x.res;
        end

        if ~isfield(x,'tunits')
            img.tunits = '';
        else
            img.tunits = x.tunits;
        end

        if ~isfield(x,'toffset')
            img.toffset = 0;
        else
            img.toffset = x.toffset;
        end

        if isfield(x,'name')
            name = x.name;
        end
                
    else
        try % Try to guess template
            [~,mask_fname] = parcellation.guess_template(x);
            [mask,img.res,img.xform,~,img.toffset,img.tunits] = nii.load(mask_fname);
            if ndims(x) < 3 % Can only reshape if we successfully guessed the template
                img.vol = matrix2vols(x,mask);
            else
                img.vol = x;
            end
            fprintf(2,'Guessing template: %s\n',mask_fname);
        catch ME
            switch ME.identifier
                case 'osl:parcellation:no_matching_mask'
                    fprintf(2,'No matching mask found - could not guess xform!\n')
                    img.xform = eye(4);
                    img.res = ones(4,1);
                    img.toffset = 0;
                    img.tunits = '';
                otherwise
                    rethrow(ME)
            end
            img.vol = x;
        end

        assert(ndims(img.vol)>=3,'Input image must be a 3D or 4D matrix, or have the same number of rows as a standard mask')
    end
    img.vol = double(img.vol);
end
