function [bad_times,final_offset] = read_bad_times_from_maxfilter(fname)
	% Parse maxfilter stdout and extract bad times
	% Output is n x 2 matrix where row is a bad segment, 
	% first column is start time, second column is stop time
	%
	% Romesh Abeysuriya 2017

	fid = fopen(fname,'r');
	l = fgetl(fid);

	initial_skip = NaN;
	extra_offset = 0;
	bad_tags = [];
	this_tag = 0;
	started_reading = false;
	reading_tag = 0; % How many tags are currently in the buffer

	while ischar(l)

		l = fgetl(fid);

		if ~isfinite(initial_skip) && ~isempty(strfind(l,'initial skip'))
			l = l(strfind(l,'#s'):end-2);
			l = strrep(l,'#s =','');
			initial_skip = str2num(l);
		end

		if strfind(l,'Reading raw tag')	
			if reading_tag == 0
				l
				this_tag = sscanf(l(strfind(l,'#t'):end),'#t = %f');
			end
			reading_tag = reading_tag + 1;
		end

		if ~isempty(strfind(l,'data block is skipped')) && ~started_reading
			extra_offset = extra_offset + 1;
		elseif ~isempty(strfind(l,'data block is skipped'))
			error('Do not know how to handle skipped blocks in the middle of the log');
		end

		if strfind(l,'Skipped') == 1
			this_tag
			l
			if ~isempty(bad_tags) && this_tag == bad_tags(end,2)
				bad_tags(end,2) = this_tag + reading_tag;
			else
				bad_tags(end+1,:) = [this_tag this_tag + reading_tag];
			end
			reading_tag = 0;
			started_reading = true;
		elseif strfind(l,'Processed raw data buffer') == 1
			reading_tag = 0;
			started_reading = true;
		end






	end

extra_offset
initial_skip

	fclose(fid);

	final_offset = - initial_skip - extra_offset; % Add this to raw Elekta times to get the final time
	bad_times = bad_tags + final_offset;
end
