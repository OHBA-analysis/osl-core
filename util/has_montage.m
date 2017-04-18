function varargout = has_montage(D,pattern,case_sensitive)
	% Check whether an MEEG has montages with names satisfying a search
	%
	% [n_present,montage_number] = has_montage(D,pattern,case_sensitive)
	%
	% INPUTS
	% - D - An MEEG object
	% - pattern - A search string (default='')
	% - case_sensitive - If true, case must match as well (default=false)
	%
	% If you run 'has_montage(D)' with no other arguments, a summary of the
	% available montages will be printed
	%
	% OUTPUTS
	% - n_present - Number of montages containing the pattern in the montage name
	% - montage_number - List of length n_present specifying the montage index of each match
	%
	% TYPICAL USAGE
	% 
	% - If you know that your analysis pipeline makes a montage with a particular
	%   name, you can check if that name is present to see if the pipeline has been run
	% - If your function requires a particular montage to be selected, you can find
	%   the required montage number using this function, and then switch to it 
	%
	% EXAMPLE
	%
	% >> has_montage(D)
	% 1 - without weights normalisation, class 1 (3559 channels)
	% 2 - with weights normalisation, class 1 (3559 channels)
	% 3 - Parcellated with weights normalisation, class 1 (68 channels)
	% >> [n,idx] = has_montage(D,'parcellated')
	%	n =
	%	     1
	%	idx =
	%	     3
	% >> D.montage('getmontage',idx).name
	%	ans =
	%	    'Parcellated with weights normalisation, class 1'

	if nargin < 3 || isempty(case_sensitive) 
		case_sensitive = false;
	end
	
	n_present = 0;
	montage_number = [];
	n_montages = D.montage('getnumber');

	if nargout > 0
		varargout{1} = n_present;
		varargout{2} = montage_number;
	end
	
	if nargin == 1
		% If no search pattern, print the montages and return
		if n_montages == 0
			fprintf('No montages present\n');
		else
			for j = 1:n_montages
				m = D.montage('getmontage',j);
				fprintf('%d - %s (%d channels)\n',j,m.name,size(m.tra,1));
			end
		end
		return
	end
	
	if n_montages == 0 % No online montages
		return
	end

	for j = 1:n_montages
		m = D.montage('getmontage',j);
		if case_sensitive
			f = strfind(m.name,pattern);
		else
			f = strfind(lower(m.name),lower(pattern));
		end

		if ~isempty(f)
			n_present = n_present + 1;
			montage_number(end+1) = j;
		end
	end

	if nargout > 0
		varargout{1} = n_present;
		varargout{2} = montage_number;
	end