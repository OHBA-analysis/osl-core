classdef (Abstract) osl_conf

	methods(Static)

		function s = read(fname)
	
			if nargin < 1 || isempty(fname) 
				fname = fullfile(osldir,'osl.conf');
			end

			if ~exist(fname,'file')
				s = struct;
				return
			end
			
			% Parse the file
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
			fclose(f);
		end

		function write(s,fname)
		
			if ~isempty(getCurrentTask)
				fprintf('Running on parallel worker, not writing to osl.conf to avoid corrupting it\n');
				return
			end

			if nargin < 2 || isempty(fname) 
				fname = fullfile(osldir,'osl.conf');
			end

			if ~isstruct(s)
				error('First input to osl_conf.write() must be a struct')
			end

			fn = fieldnames(s);
			f = fopen(fname,'w');

			for j = 1:length(fn)
				fprintf(f,'%s=%s\n',fn{j},s.(fn{j}));
			end

			fclose(f);

		end

	end

end


