function s = read_conf_file(fname)
	% Read ASCII config file with key-value pairs
	% Each line consists of 3 parts
	% - key: A valid MATLAB struct field name e.g. doesn't start with a number
	% - equals sign (separating key and value)
	% - value: String that is assigned to struct
	% Note that leading and trailing spaces in the key and value will be stirpped
	%
	% RA 2017
	
	if nargin < 1 || isempty(fname) 
		fname = fullfile(osldir,'osl.conf');
	end

	f = fopen(fname,'r');
	l = fgetl(f);
	s = struct();
	while l ~= -1
        if ~isempty(strtrim(l))
            lp = regexp(l,'=','split');
            if length(lp)>2
            	lp{2} = sprintf('%s=',lp{2:end});
            	lp{2} = lp{2}(1:end-1);
            end
            s.(strtrim(lp{1})) = strtrim(lp{2});
        end
        l = fgetl(f);    
	end




	

