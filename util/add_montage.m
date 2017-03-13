function D = add_montage(D,W,name,ROIlabels)
	% Add an online montage to an MEEG object (in memory)
	% 
	% INPUTS
	% - D (MEEG object)
	% - W (transformation matrix being applied to current montage (n_new_channels x n_old_channels)
	% 	  This matrix premultiplies the active montage and must therefore have the same number 
	
	if nargin < 4 || isempty(ROIlabels) 
		ROIlabels = [];
	end
	
	if nargin < 3 || isempty(name) 
		name = 'New montage';
	end
	
	nMontages       = D.montage('getnumber');
	currentMontage  = D.montage('getmontage');

	assert(size(W,2)==size(currentMontage.tra,1),sprintf('Transformation matrix is wrong size (is %dx%d, needs to be Nx%d)',size(W,1),size(W,2),size(currentMontage.tra(1))))

	newMontage      = currentMontage;
	unit            = unique(D.units());
	newMontage.name = name;

	if isempty(newMontage) 
		% If no online montage is currently applied, use weights directly
		newMontage.tra  = W;
	else
		newMontage.tra  = W * currentMontage.tra;
	end

	if ~isempty(ROIlabels)
		newMontage.labelnew = ROIlabels;
	end

	D = D.montage('add', newMontage);
	D = D.montage('switch', nMontages + 1);
	D = D.units(:,unit{1});
