function bad_times = read_bad_times_from_maxfilter(fname)
	% Parse maxfilter stdout and extract bad times
	% Output is n x 2 matrix where row is a bad segment, 
	% first column is start time, second column is stop time
	%
	% Romesh Abeysuriya 2017

	fid = fopen(fname,'r');
	l = fgetl(fid);

	initial_skip = NaN;
	computed_offset = NaN;
	new_block = true; % This flag lets us perform operations once per block
	bad_tags = [];

	while ischar(l)

		l = fgetl(fid);

		if ~isfinite(initial_skip) && ~isempty(strfind(l,'initial skip'))
			l = l(strfind(l,'#s'):end-2);
			l = strrep(l,'#s =','');
			initial_skip = str2num(l);
		end

		if strfind(l,'Reading raw tag')	
			this_tag = sscanf(l(strfind(l,'#t'):end),'#t = %f');
			new_block = true;
		end

		if new_block & strfind(l,'Time ') == 1
			this_offset = this_tag - sscanf(l,'Time %f');

			if ~isfinite(computed_offset)
				computed_offset = this_offset;
			else
				assert(this_offset == computed_offset);
			end

			new_block = false;
		elseif new_block & strfind(l,'#t = ') == 1
			this_offset = this_tag - sscanf(l,'#t = %f');

			if ~isfinite(computed_offset)
				computed_offset = this_offset;
			else
				assert(this_offset == computed_offset);
			end

			new_block = false;
		end

		% --- Reading raw tag #b = 3/1311 (#t = 83.000) ---

		% Time 81.000: cont HPI is off, data block is skipped!



		if strfind(l,'Skipped 1 bad data tags')
			if ~isempty(bad_tags) && this_tag == bad_tags(end,2)
				bad_tags(end,2) = this_tag + 1;
			else
				bad_tags(end+1,:) = [this_tag this_tag + 1];
			end
		end

	end

	fclose(fid);

	bad_times = bad_tags - initial_skip - computed_offset;
end
