function D = add_montage(D,W,name,labels,channels)
	% Add an online montage to an MEEG object (in memory)
	% 
	% INPUTS
	% - D (MEEG object)
	% - W (transformation matrix being applied to current montage (n_new_channels x n_old_channels)
	% 	  This matrix premultiplies the active montage and must therefore have the same number 
	% - name (optional string giving name of new montage)
	% - labels (optional cell array of names for each channel)
	% - channels (optional struct array that provides channel information)

	% First, check proposed tra matrix is OK in terms of orientation
	% That way, size(W,1) is guaranteed to be the new number of channels
	assert(size(W,2)==D.nchannels,sprintf('Transformation matrix is wrong size (is %dx%d, needs to be Nx%d)',size(W,1),size(W,2),D.nchannels));

	if nargin < 5 || isempty(channels) 
		channels = [];
	end
	
    if nargin < 4 || isempty(labels) 
        labels = arrayfun(@(x) sprintf('ROI%d',x),1:size(W,1),'UniformOutput',false);
    else
    	assert(iscell(labels) && length(labels)==size(W,1),sprintf('Channel labels must be a cell array with %d elements',size(W,1)));
	end

	if nargin < 3 || isempty(name) 
		name = 'New montage';
	end

	currentMontage  = D.montage('getmontage');
    
    if isempty(currentMontage)
    	newMontage = struct; 
    	newMontage.tra  = W;
    	newMontage.labelorg = D.chanlabels;
    else
    	newMontage      = currentMontage;
		newMontage.tra  = W * currentMontage.tra;
	end

	if isempty(channels)
        if isfield(newMontage,'channels')
            newMontage      = rmfield(newMontage, 'channels');
        end
	else
		newMontage.channels = channels;
	end
	
	newMontage.name = name;
	newMontage.labelnew = labels;

	D = D.montage('add', newMontage);
	D = D.montage('switch', D.montage('getnumber'));
	