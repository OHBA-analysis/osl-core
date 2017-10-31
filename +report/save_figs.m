function save_figs(h,dirname)
	% Take in an array of figure handles
	% Save them to directory 'dirname'
	% Create directory if required
	% Automatically infer filename if required

	if ~isdir(dirname)
		mkdir(dirname);
	end

	for j = 1:length(h)
		fname = get(h(j),'tag');
		if isempty(fname)
			[~,tmp] = fileparts(tempname('.'));
			fname = ['untitled_' tmp(1:10)];
		end
		
		pfig(h(j),fullfile(dirname,[fname '.png']));
	end
