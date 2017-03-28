function D = add_montage(D,W,name,labels)
	% Add an online montage to an MEEG object (in memory)
	% 
	% INPUTS
	% - D (MEEG object)
	% - W (transformation matrix being applied to current montage (n_new_channels x n_old_channels)
	% 	  This matrix premultiplies the active montage and must therefore have the same number 
	% - name (optional string giving name of new montage)
	% - labels (optional cell array of names for each channel)
	
	if nargin < 4 || isempty(labels) 
		labels = [];
	end
	
	if nargin < 3 || isempty(name) 
		name = 'New montage';
	end
	
	nMontages       = D.montage('getnumber');
	currentMontage  = D.montage('getmontage');

	if size(W,2) ~= size(currentMontage.tra,1) && size(W,1) == size(currentMontage.tra,1)
		W = W.';
		fprintf(2,'Warning - input transformation matrix is being transposed. It was %d x %d, and the MEEG montage has %d channels. To avoid this message, transpose the input to add_montage()\n',size(W,1),size(W,2),size(currentMontage.tra,1));
    end
    
    if nargin < 4 || isempty(labels) 
        labels = arrayfun(@(x) sprintf('ROI%d',x),1:size(W,1),'UniformOutput',false);
	end

	assert(size(W,2)==size(currentMontage.tra,1),sprintf('Transformation matrix is wrong size (is %dx%d, needs to be Nx%d)',size(W,1),size(W,2),size(currentMontage.tra,1)))
	assert(iscell(labels) && length(labels)==size(W,1),sprintf('ROI labels must be a cell array with %d elements',size(W,1)));

	newMontage      = currentMontage;
    newMontage      = rmfield(newMontage, 'channels');
	unit            = unique(D.units());
	newMontage.name = name;
    
	if isempty(newMontage) 
		% If no online montage is currently applied, use weights directly
		newMontage.tra  = W;
	else
		newMontage.tra  = W * currentMontage.tra;
	end

	if ~isempty(labels)
		newMontage.labelnew = labels;
    end

	D = D.montage('add', newMontage);
	D = D.montage('switch', nMontages + 1);
	D = D.units(:,unit{1});
