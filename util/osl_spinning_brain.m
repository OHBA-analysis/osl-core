function osl_spinning_brain(fname,fig,duration)
	% Make an animated gif of a spinning brain
	% The azimuth angle will sweep through 1 revolution, while the
	% elevation will remain unchanged
	% The animation will start with the current azimuth angle
	%
	% INPUTS
	% - fname - Input file
	% - fig - Handle to figure to rotate. All axes in figure will be rotated
	% - duration - the period of 1 rotation in seconds
	% - elevation - Override the elevation in the plots

	if nargin < 3 || isempty(duration) 
		duration = 5;
	end
	
	if nargin < 2 || isempty(fig) 
		fig = gcf;
	end

	if nargin < 1 || isempty(fname) 
		error('Must specify an output filename')
	end
	

	fps = 25;
	delayTime = 1/fps;
	n_frames = duration*fps;

	ax = findobj(fig,'Type','axes'); 
	axis(ax,'vis3d')

	for j = 1:n_frames

		for k = 1:length(ax)
			camorbit(ax(k),360/n_frames,0,'data',[0 0 1])
		end

		cdata = getframe(fig);
		im = frame2im(cdata);

		%cdata = print('-RGBImage');
		[A, map] = rgb2ind(im, 256, 'nodither');
		if j == 1
			imwrite(A, map, fname, 'gif', 'DelayTime', delayTime, 'LoopCount', Inf);
		else
			imwrite(A, map, fname, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
		end
	end



