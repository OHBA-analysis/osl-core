classdef (Abstract) osl_conf

	methods(Static)

		function s = read(fname)
	
			if nargin < 1 || isempty(fname) 
				fname = fullfile(osldir,'osl.conf');
			end

			if ~isfile(fname)
				s = struct();
				return
			end
			
			% Parse the file
			f = fopen(fname,'r');
			l = fgetl(f);
			s = struct();
			while l ~= -1
		        if ~isempty(strtrim(l))
		            [name,value] = strtok(l,'=');
		            s.(strtrim(name)) = strtrim(value);
		        end
		        l = fgetl(f);    
			end
			fclose(f);

		end

		function write(t,fname)
		
			if ~isempty(getCurrentTask)
				fprintf('Running on parallel worker, not writing to osl.conf to avoid corrupting it\n');
				return
			end

			if nargin < 2 || isempty(fname) 
				fname = fullfile(osldir,'osl.conf');
			end

			assert( isstruct(t) && isscalar(t), 'First input to osl_conf.write() must be a struct.' );
			assert( all(structfun( @ischar, t )), 'Config values must be strings.' );
			
			% load existing config and overwrite fields
			s = osl_conf.read(fname);
			s = structmerge( s, t );
			
			f = fopen(fname,'w+');
			structkvfun( @(k,v) fprintf(f,'%s = %s\n',k,v), s );
			fclose(f);

		end

	end

end