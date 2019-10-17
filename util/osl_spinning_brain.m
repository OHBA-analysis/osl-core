function osl_spinning_brain(fname,fig,duration,tight)
    % Make an animation by spinning the figure
    % The azimuth angle will sweep through 1 revolution, while the
    % elevation will remain unchanged
    % The animation will start with the current azimuth angle
    % All axes will be rotated
    %
    % INPUTS
    % - fname - output file name. Format is chosen based on extension - can be '.mp4' or '.gif'
    % - fig - Handle to figure to rotate. All axes in figure will be rotated (default: current figure)
    % - duration - the period of 1 rotation in seconds (default: 5 seconds)
    % - tight - Enlarge brain in plot by expanding the axis and reducing whitespace. Colorbar will disappear!

    % FILE FORMAT NOTE
    %
    % mp4 videos are much less CPU-intensive than gif animations, particularly when there are
    % multiple videos on a single slide. To animate:
    %
    % KEYNOTE - Uncheck 'Start movie on click' and set 'Repeat' to 'Loop'
    % POWERPOINT - Check 'Loop until stopped' and add an animation with emphasis effect 'Play' to start
    % BROWSER HTML - "<video src="filename.mp4" autoplay loop>"

    if nargin < 4 || isempty(tight) 
        tight = false;
    end
    
    if nargin < 3 || isempty(duration) 
        duration = 5;
    end
    
    if nargin < 2 || isempty(fig) 
        fig = gcf;
    end

    if nargin < 1 || isempty(fname) 
        error('Must specify an output filename')
    end
    
    [~,~,ext] = fileparts(fname);
    if ~any(strcmp(ext,{'.gif','.mp4'}))
        error('File name should end in .gif or .mp4');
    end

    fps = 30;
    n_frames = duration*fps;

    ax = findobj(fig,'Type','axes'); 
    axis(ax,'vis3d')

    if tight
        set(ax,'Position',[0 0 1 1]);
        outerpos = ax.OuterPosition;
        ax.Position = [outerpos(1) outerpos(2) 1-2*outerpos(1) 1-2*outerpos(2)];
    end

    if strcmp(ext,'.mp4')
        v = VideoWriter(fname,'MPEG-4');
        v.FrameRate = fps;
        v.Quality = 60; % About half the default quality file size with minimal loss of quality
        open(v)
        cleanup = onCleanup(@() close(v)); % Automatically close when function ends
    else
        delayTime = 1/fps;
    end

    for j = 1:n_frames
        
        for k = 1:length(ax)
            current_view = get(ax(k),'View');
            camorbit(ax(k),sign(current_view(2)+eps)*360/n_frames,0);
        end

        cdata = getframe(fig);
        im = frame2im(cdata);

        if strcmp(ext,'.mp4')
            writeVideo(v,im);
        else
            [A, map] = rgb2ind(im, 256, 'nodither');
            if j == 1
                imwrite(A, map, fname, 'gif', 'DelayTime', delayTime, 'LoopCount', Inf);
            else
                imwrite(A, map, fname, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
            end
        end
    end


