classdef osleyes < handle
	% A pure Matlab NIFTI viewer based on common OSL usage of fsleyes/fslview
	%
	% 
	% TERMINOLOGY
	% - A *layer* corresponds to a NIFTI file
	% - A *volume* corresponds to a slice in the 4th dimension of the data in the NIFTI file
	%
	% osleyes.colormaps - A cell array of colormaps for each layer. Each
	% 					  element corresponds to a layer. The layer can
	% 					  contain a colormap, or a cell array with 2 colormaps
	% 					  for a bidirectional map. A colormap is 
	% 						- the name of a function that takes in a number of elements e.g. 'hot'
	%						- a (N x 3) matrix of RGB values (e.g. the matrix returned by 'osl_colormap('hot')'
	%
	% OTHER NOTES
	% - figure contains 'osleyes' property with handle to this object
	% - figure will be closed if handle is deleted
	% - object will be deleted if figure is closed

	% Romesh Abeysuriya 2017

	properties
		clims = {}; % Values below this range are transparent, values above are clipped
		colormaps= {}; % Can be name of a colormap function on the path
		current_point = [1 1 1]; % These are the XYZ MNI coordinates for the crosshair
		current_vols = []; % This is an array that specifies which volume is being displayed for each layer
		visible = logical(1); % This is an array that is true/false for whether a layer is visible or not
		active_layer = 1; % The active layer controls which layer's properties are shown in the control panel and the colorbars
		show_controls = 1;
		show_crosshair = 1;
		layer_alpha = []; % Transparency of each layer
	end

	properties(SetAccess=protected)
		niifiles
		img = {} % Temporarily accessible
		xform
		colormap_resolution = 255; % Number of colors in each colormap (if not specified as matrix)
	end

	properties(GetAccess=private,SetAccess=private)
		fig % Handle to the main figure window bound to the object
		ax % Handles of the three display axes
		controls % Handles for control panel and associated controls

		ts_ax  % Handle to timeseries axis
		ts_line  % Handle to line in timeseries plot
		ts_bar  % Handle to marker bar in timeseries plot

		h_img = {} % Cell array of img handles associated with each nii file (there are 3)
		h_crosshair % Handles for crosshairs on each axis
		h_coloraxes % Handles to colorbar axes
		h_colorimage % Handles to colorbar images
		
		coord % Axis coordinates for each image
		colormap_matrices = {}; % Cached colormaps
		under_construction = true; % Don't render anything while this is true
		motion_active = 0; % If this flag is true, then the slices will be updated when the mouse is moved
		lims = nan(3,2); % Axis limits for each MNI dimension (used in plots)
		contextmenu % Handle for the context menu
	end

	methods

		function self = osleyes(niifiles,clims,colormaps)
			% niifiles is a file name or cell array of file names
			% If it's a cell array with an empty first element, then the mask
			% will be automatically guessed

			if nargin < 1 || isempty(niifiles) 
				niifiles = {fullfile(osldir,'std_masks/MNI152_T1_2mm_brain.nii.gz')};
            end
			
            if ~iscell(niifiles)
                niifiles = {[],niifiles};
            end
            
            if nargin < 3 || isempty(colormaps) 
            	colormaps = cell(length(niifiles),1);
            end

            if nargin < 2 || isempty(clims) 
            	clims = cell(length(niifiles),1);
            end

            if isempty(niifiles{1});
            	vol = nii.load(niifiles{2});
            	try
            		[~,niifiles{1},~] = parcellation.guess_template(vol);
            		clims = {[] clims{:}};
            		colormaps = {[] colormaps{:}};
            	catch ME
            		ME.getReport
            		niifiles = niifiles(2:end);
            		clims = clims(2:end);
            		colormaps = colormaps(2:end);
            	end
        	end

			self.fig = figure('Units','Characters','Color','k');
			self.initial_render();
			set(self.fig,'KeyPressFcn',@(a,b) KeyPressFcn(self,a,b),'CloseRequestFcn',@(~,~) delete(self),'ResizeFcn',@(~,~) resize(self));
			addprop(self.fig,'osleyes');
			set(self.fig,'osleyes',self); % Store handle to this osleyes in the figure so it can be retrieved later if desired

			self.niifiles = niifiles;
			dropdown_strings = {};
			for j = 1:length(self.niifiles)
				[self.img{j},~,self.xform{j}] = nii.load(self.niifiles{j}); % Do not apply xform/qform
				self.img{j} = double(self.img{j});
				[~,fname,ext] = fileparts(self.niifiles{j});
				dropdown_strings{j} = [fname '.' ext];
				self.coord{j} = get_coords(self.img{j},self.xform{j});
				self.h_img{j}(1) = image(self.coord{j}.y,self.coord{j}.z,permute(self.img{j}(1,:,:,1),[2 3 1])','Parent',self.ax(1),'HitTest','off','AlphaDataMapping','none','AlphaData',1);
				self.h_img{j}(2) = image(self.coord{j}.x,self.coord{j}.z,permute(self.img{j}(:,1,:,1),[1 3 2])','Parent',self.ax(2),'HitTest','off','AlphaDataMapping','none','AlphaData',1);
				self.h_img{j}(3) = image(self.coord{j}.x,self.coord{j}.y,permute(self.img{j}(:,:,1,1),[1 2 3])','Parent',self.ax(3),'HitTest','off','AlphaDataMapping','none','AlphaData',1);
				
				if isempty(clims{j})
					self.clims{j} = [0,max(self.img{j}(:))];
				else
					self.clims{j} = clims{j};
				end

				if isempty(colormaps{j})
					if min(self.img{j}(:)) < 0
						self.colormaps{j} = {osl_colormap('hot'),osl_colormap('cold')};
					elseif j == 1
						self.colormaps{j} = osl_colormap('grey');
					else
						self.colormaps{j} = osl_colormap('hot');
					end
				else
					self.colormaps{j} = colormaps{j};
				end

				self.current_vols(j) = 1;
				self.visible(j) = 1;
				self.layer_alpha(j) = 1;
			end

			set(self.controls.image_list,'String',dropdown_strings);
			self.lims(1,:) = [min(cellfun(@(x) min(x.x),self.coord)) max(cellfun(@(x) max(x.x),self.coord))];
			self.lims(2,:) = [min(cellfun(@(x) min(x.y),self.coord)) max(cellfun(@(x) max(x.y),self.coord))];
			self.lims(3,:) = [min(cellfun(@(x) min(x.z),self.coord)) max(cellfun(@(x) max(x.z),self.coord))];
			
			set(self.ax(1),'XLim',self.lims(2,:),'YLim',self.lims(3,:),'Visible','off');
			set(self.ax(2),'XLim',self.lims(1,:),'YLim',self.lims(3,:),'Visible','off');
			set(self.ax(3),'XLim',self.lims(1,:),'YLim',self.lims(2,:),'Visible','off');

			set(self.ax(1),'Color','k','Clipping','off','View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal');
			set(self.ax(2),'Color','k','Clipping','off','View',[0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
			set(self.ax(3),'Color','k','Clipping','off','View', [0 90],'DataAspectRatio',[1 1 1],'YDir','normal','XDir','reverse');
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
			self.refresh_colors();
			self.active_layer = length(self.niifiles);
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
			
			if size(self.img{self.active_layer},4)==1
				error('Layer only has one volume, there is nothing to animate!');
			end

			while 1
				self.current_vols(self.active_layer) = 1+mod(self.current_vols(self.active_layer),size(self.img{self.active_layer},4));
				pause(1/fps);
			end
		end
		
		function set.colormaps(self,val)
			if self.under_construction
				self.colormaps = val;
				return;
			end

			assert(iscell(val),'Colormaps should be a cell array')
			assert(length(val) == length(self.img),'Number of colormaps must match number of images');
			for j = 1:length(val)
				if iscell(val{j})
					assert(length(val{j})==2,'Colormap must either be a value or a cell array of length 2');
				end
			end
			self.colormaps = val;
			self.refresh_colors;
		end

		function set.clims(self,val)
			if self.under_construction
				self.clims = val;
				return;
			end

			assert(iscell(val),'clims should be a cell array')
			assert(length(val) == length(self.img),'Number of clims must match number of images');
			for j = 1:length(self.colormaps)
				assert(length(val{j})==2,'Colour limits must have two elements')
				assert(val{j}(2)>=val{j}(1),'Colour limits must be in ascending order')
				if iscell(self.colormaps{j})
					assert(val{j}(1)>=0,'If using bidirectional limits, colour range must start >= 0')
				end
			end
			self.clims = val;
			self.refresh_colors;
			self.active_layer = self.active_layer; % Update the text boxes and colorbar limits
		end

		function set.current_point(self,val)
			% Validate limits
			xl = get(self.ax(3),'XLim');
			yl = get(self.ax(1),'XLim');
			zl = get(self.ax(1),'YLim');

			if val(1) > xl(2) || val(1) < xl(1)
				return
			end

			if val(2) > yl(2) || val(2) < yl(1)
				return
			end

			if val(3) > zl(2) || val(3) < zl(1)
				return
			end

			self.current_point = val;
			refresh_slices(self)
		end

		function v = set.active_layer(self,val)
			assert(val > 0,'Layer must be positive');
			assert(val <= length(self.niifiles),'Layer number cannot exceed number of layers');
			assert(val == round(val),'Layer must be an integer');

			self.active_layer = val;
			set(self.controls.image_list,'Value',val)
			set(self.controls.clim(1),'String',sprintf('%g',self.clims{get(self.controls.image_list,'Value')}(1)));
			set(self.controls.clim(2),'String',sprintf('%g',self.clims{get(self.controls.image_list,'Value')}(2)));
			set(self.controls.volume,'String',sprintf('%d',self.current_vols(get(self.controls.image_list,'Value'))));
			set(self.controls.volume_label_count,'String',sprintf('of %d',size(self.img{self.active_layer},4)));

			set(self.controls.visible,'Value',self.visible(self.active_layer));

			tickstrs = @(low,high,n) arrayfun(@(x) sprintf('%.2g',x),linspace(low,high,n),'UniformOutput',false);

			if iscell(self.colormaps{self.active_layer})
				set(self.h_coloraxes(1),'Visible','on');
				set(self.h_colorimage(1),'Visible','on');
				set(self.h_colorimage(2),'CData',permute(self.colormap_matrices{self.active_layer}{1},[1 3 2]));
				set(self.h_colorimage(1),'CData',permute(self.colormap_matrices{self.active_layer}{2},[1 3 2]));
				set(self.h_coloraxes(1),'YTick',linspace(0,1,4),'YTickLabel',tickstrs(-self.clims{self.active_layer}(1),-self.clims{self.active_layer}(2),4))
				set(self.h_coloraxes(2),'YTick',linspace(0,1,4),'YTickLabel',tickstrs(self.clims{self.active_layer}(1),self.clims{self.active_layer}(2),4))
			else
				set(self.h_coloraxes(1),'Visible','off');
				set(self.h_colorimage(1),'Visible','off');
				set(self.h_colorimage(2),'CData',permute(self.colormap_matrices{self.active_layer},[1 3 2]));
				set(self.h_coloraxes(2),'YTick',linspace(0,1,4),'YTickLabel',tickstrs(self.clims{self.active_layer}(1),self.clims{self.active_layer}(2),4))
			end

			if size(self.img{self.active_layer},4)==1
				set(self.controls.volume,'Enable','off');
			else
				set(self.controls.volume,'Enable','on');
			end

			self.resize(); % Update sizes of colorbars
			self.refresh_slices();
		end

		function set.current_vols(self,val)
			if self.under_construction
				self.current_vols = val;
				return;
			end

			assert(isvector(val) && length(val) == length(self.niifiles));

			for j = 1:length(val)
				assert(val(j)>0 && val(j) <= size(self.img{j},4),'Volume index out of bounds');
			end

			self.current_vols = val;
			set(self.controls.volume,'String',sprintf('%d',self.current_vols(get(self.controls.image_list,'Value'))));
			self.refresh_slices();

			if ishandle(self.ts_bar)
				set(self.ts_bar,'XData',[1 1]*self.current_vols(self.active_layer),'YData',get(get(self.ts_bar,'Parent'),'YLim'));
			end

		end

		function set.visible(self,val)
			self.visible = val;

			if self.under_construction
				return;
			end

			self.refresh_slices();
		end

		function set.show_controls(self,val)
			self.show_controls = logical(val);
			self.resize()
		end

		function set.show_crosshair(self,val)
			self.show_crosshair = logical(val);
			self.refresh_slices()
		end

		function set.layer_alpha(self,val)
			self.layer_alpha = val;

			if self.under_construction
				return;
			end
			
			self.refresh_slices()
		end

	end

	methods(Access=private)

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
			ax = axes(fig);
			self.ts_line = plot(ax,1,1,'HitTest','off');
			hold(ax,'on')
			self.ts_bar = plot(ax,1,1,'r','HitTest','off');
			set(ax,'ButtonDownFcn',@(~,~) set_volume(ax,self));
			self.current_vols = self.current_vols; % Reset the data in h_bar via set.current_vols

		end

		function refresh_slices(self)
			% Update the slices
			% This will move the crosshairs, set the position markers, and change which slice is displayed
			% It will not re-render the colour bars

			p = self.current_point; % Current point in 3D

			if self.show_crosshair
				set(self.h_crosshair(1),'Visible','on','XData',[self.lims(2,:) NaN p(2) p(2)],'YData',[p(3) p(3) NaN self.lims(3,:)]);
				set(self.h_crosshair(2),'Visible','on','XData',[self.lims(1,:) NaN p(1) p(1)],'YData',[p(3) p(3) NaN self.lims(3,:)]);
				set(self.h_crosshair(3),'Visible','on','XData',[self.lims(1,:) NaN p(1) p(1)],'YData',[p(2) p(2) NaN self.lims(2,:)]);
			else
				set(self.h_crosshair(1),'Visible','off');
				set(self.h_crosshair(2),'Visible','off');
				set(self.h_crosshair(3),'Visible','off');
			end

			set(self.controls.marker(1),'String',sprintf('X = %+06.1f',p(1)));
			set(self.controls.marker(2),'String',sprintf('Y = %+06.1f',p(2)));
			set(self.controls.marker(3),'String',sprintf('Z = %+06.1f',p(3)));
			%self.controls.marker(4)

			% Now update each slice

			for j = 1:length(self.img)
				[~,idx(1)] = min(abs(self.coord{j}.x-p(1)));
				[~,idx(2)] = min(abs(self.coord{j}.y-p(2)));
				[~,idx(3)] = min(abs(self.coord{j}.z-p(3)));

				% These are the slice maps - need to convert them to color values now
				d1 = permute(self.img{j}(idx(1),:,:,self.current_vols(j)),[2 3 1])';
				d2 = permute(self.img{j}(:,idx(2),:,self.current_vols(j)),[1 3 2])';
				d3 = permute(self.img{j}(:,:,idx(3),self.current_vols(j)),[1 2 3])';

				if j == self.active_layer
					set(self.controls.marker(4),'String',sprintf('Value = %+ 7.2f',self.img{j}(idx(1),idx(2),idx(3),self.current_vols(j))));
                end


                % If bidirectional colormap, then hide abs(x)<clim(1), otherwise hide x<clim(1)
                if iscell(self.colormaps{j})
                	hidefcn = @abs; 
                else
                	hidefcn = @(x) x;
                end

                if self.visible(j)
					set(self.h_img{j}(1),'Visible','on','CData',map_colors(d1,self.colormap_matrices{j},self.clims{j}),'AlphaData',self.layer_alpha(j)*+(hidefcn(d1)>self.clims{j}(1)));
					set(self.h_img{j}(2),'Visible','on','CData',map_colors(d2,self.colormap_matrices{j},self.clims{j}),'AlphaData',self.layer_alpha(j)*+(hidefcn(d2)>self.clims{j}(1)));
					set(self.h_img{j}(3),'Visible','on','CData',map_colors(d3,self.colormap_matrices{j},self.clims{j}),'AlphaData',self.layer_alpha(j)*+(hidefcn(d3)>self.clims{j}(1)));
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
				set(self.ts_line,'XData',1:size(self.img{self.active_layer},4),'YData',squeeze(self.img{self.active_layer}(idx(1),idx(2),idx(3),:)));
				title(get(self.ts_line,'Parent'),sprintf('%s - MNI (%.2f,%.2f,%.2f)',self.controls.image_list.String{self.active_layer},p(1),p(2),p(3)),'Interpreter','none')
			end
		end

		function refresh_colors(self)
			% Turn the colormap strings into colormap matrices
			% Called if clims or colormap is changed
			for j = 1:length(self.colormaps)
				self.colormap_matrices{j} = self.compute_color_matrix(self.colormaps{j});
			end

			% Setting the value here causes the image to be refreshed with the new colour maps
			% active_layer is changed because the colour bars also need to be refreshed
			self.active_layer = self.active_layer;
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

			%self.h_colorbar(1) = colorbar('Peer',self.ax(3),'Location','eastoutside','Color','w','Units','characters','TickDirection','in','AxisLocation','in');
			%self.h_colorbar(2) = colorbar('Peer',self.ax(3),'Location','westoutside','Color','w','Units','characters','TickDirection','in','AxisLocation','out');
			self.h_coloraxes(1) = axes('Box','on','Color','k','Units','characters');
			self.h_coloraxes(2) = axes('Box','on','Color','k','Units','characters');
			self.h_colorimage(1) = image([0 1],[0 1],1,'Parent',self.h_coloraxes(1));
			self.h_colorimage(2) = image([0 1],[0 1],1,'Parent',self.h_coloraxes(2));
			set(self.h_coloraxes,'XLim',[0 1],'YLim',[0 1],'XColor','w','YColor','w','XTick',[],'YDir','reverse','YAxisLocation','right');
			set(self.h_coloraxes(2),'YDir','normal');
			
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
			self.controls.marker(4) = uicontrol(self.controls.panel,'style','text','String','Value = +0000.0','Units','characters','Position',[0 0.9 18 1],'HorizontalAlignment','left');

			% Put b to the right of a with given padding
			next_to = @(b,a,padding) 	set(b,'Position',sum(get(a,'Position').*[1 0 1 0]).*[1 0 0 0]  + [padding 0 0 0] + [0 1 1 1].*get(b,'Position'));

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

		end


		function resize(self)
			% Resize the main objects in the figure window
			% The contents of the control panel will not be changed - their position is determined
			% by the initial layout function and is fixed from that point on

			figpos = get(self.fig,'Position');

			% Layout is - 0.3 reserved for each axis, plus 0.1 for the two colorbars
			cb = 0.1; % Colorbar width			
			w = 0.3*(1-cb)*figpos(3);
			p = ((1-cb)*figpos(3)-3*w)/6;
			control_height = 3*+self.show_controls; % Control panel height in characters
			ax_height = figpos(4)-control_height;

			if ax_height <= 0
				return
			end

			set(self.ax(1),'Position',[p control_height w ax_height]);
			set(self.ax(2),'Position',[p+w+p+p  control_height w ax_height]);
			set(self.ax(3),'Position',[p+w+p+p+w+p+p control_height w ax_height]);
			set(self.controls.panel,'Position',[0 0 figpos(3) control_height]);

			if iscell(self.colormaps{self.active_layer})
				set(self.h_coloraxes(1),'Position',[(1-cb+0.005)*figpos(3) control_height+ax_height*0.075 0.015*figpos(3) ax_height*0.4]);
				set(self.h_coloraxes(2),'Position',[(1-cb+0.005)*figpos(3) control_height+ax_height*0.525 0.015*figpos(3) ax_height*0.4]);
			else
				set(self.h_coloraxes(2),'Position',[(1-cb+0.005)*figpos(3) control_height+ax_height*0.1 0.015*figpos(3) ax_height*0.8]);
			end

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
	self.clims{self.active_layer} = [str2double(get(self.controls.clim(1),'String')) str2double(get(self.controls.clim(2),'String'))];
end

function visible_box_callback(self)
	% This callback is run if the user clicks the visibility checkbox
	self.visible(self.active_layer) = get(self.controls.visible,'Value');
end

function volume_box_callback(self)
	% If the user types a umber into the volume edit text box, this callback is run
	new_volume = str2double(get(self.controls.volume,'String'));
	if new_volume < 1 || new_volume > size(self.img{self.active_layer},4)
		uiwait(errordlg(sprintf('Volume out of bounds (image has %d volumes)',size(self.img{self.active_layer},4))));
		set(self.controls.volume,'String',num2str(self.current_vols(self.active_layer)));
	else
		self.current_vols(self.active_layer) = new_volume;
	end
end

function rgb = map_colors(x,cmap,clim)
	% Map the image to colours using the colour matrices
	% The problem is that each axis contains multiple images each with 
	% a different colourmap. But an axis can only have one colormap. Therefore
	% we need to perform the colour mapping ourselves
	% 
	% INPUT
	% 	x - image data
	% 	cmap - color matrix OR cell with two color matrices
	% 	clim - color limit for this plot
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
		r(x>=0) = cmap{1}(pos_idx(x>=0),1);
		g(x>=0) = cmap{1}(pos_idx(x>=0),2);
		b(x>=0) = cmap{1}(pos_idx(x>=0),3);
		r(x<0) = cmap{2}(neg_idx(x<0),1);
		g(x<0) = cmap{2}(neg_idx(x<0),2);
		b(x<0) = cmap{2}(neg_idx(x<0),3);
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

function set_volume(ax,h_osleyes)
	% This is the callback for the axis object in the timeseries plot
	% It is wrapped in an anonymous function closure that contains
	% - ax (axes of the timeseries plot)
	% - h_bar (vertical red line showing which time was clicked)
	% - h_osleyes (handle to osleyes object that created and is thus bound to this plot)
	p = get(ax,'CurrentPoint');
	h_osleyes.current_vols(h_osleyes.active_layer) = round(p(1,1));
end

function KeyPressFcn(self,~,KeyData)
	switch KeyData.Key
		case 'uparrow'
			self.current_vols(self.active_layer) = 1+mod(self.current_vols(self.active_layer),size(self.img{self.active_layer},4));
		case 'downarrow'
			self.current_vols(self.active_layer) = 1+mod(self.current_vols(self.active_layer)-2,size(self.img{self.active_layer},4));
	end
end


