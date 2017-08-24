function D2 = osl_filter(D,fband,order)
	% Apply temporal filter, return copy

	if nargin < 3 || isempty(order) 
		order = [];
	end
	
	S = struct;

	montage_idx = D.montage('getindex');
	S.D = D.montage('switch',0);

	if fband(1) == 0
		S.band = 'low';
		S.freq = fband(2);
	elseif fband(2) == inf
		S.band = 'high';
		S.freq = fband(1);
	else
		S.band = 'bandpass';
		S.freq = fband;
	end

	if ~isempty(order)
		S.order = order;
	end

	D2 = spm_eeg_filter(S);
	D2 = D2.montage('switch',montage_idx);
	



