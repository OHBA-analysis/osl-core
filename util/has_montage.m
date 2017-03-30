function [n_present,montage_number] = has_montage(D,pattern,case_sensitive)
	% Check whether an MEEG has montages with names satisfying a search
	%
	% INPUTS
	% - D - An MEEG object
	% - pattern - A search string
	% - case_sensitive - If true, case must match as well (default=false)
	%
	% OUTPUTS
	% - n_present - Number of montages containing the pattern in the montage name
	% - montage_number - List of length n_present specifying the montage number of the matches
	%
	% EXAMPLE
	% [n,idx] = has_montage(D,'parcellated')
	%	n =
	%	     1
	%	idx =
	%	     3
	% 
	% D.montage('getmontage',idx).name
	%	ans =
	%	    'Parcellated with weights normalisation, class 1'
	
	if nargin < 3 || isempty(case_sensitive) 
		case_sensitive = false;
	end
	
	n_present = 0;
	montage_number = [];

	n_montages = D.montage('getnumber');
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
