classdef (Abstract) osl_conf

	methods(Static)

		function s = read()
	
			fname = getenv('OSLCONF');
			assert( ~isempty(fname), 'Missing path to config file (env:OSLCONF).' );

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
					parts = strsplit(l,'=');
					name  = strtrim(parts{1});
					value = strjoin(parts(2:end),'=');
		            s.(name) = strtrim(value);
		        end
		        l = fgetl(f);    
			end
			fclose(f);

		end

		function write(s)
		
			if ~isempty(getCurrentTask)
				fprintf('Running on parallel worker, not writing to osl.conf to avoid corrupting it\n');
				return
			end

			fname = getenv('OSLCONF');
			assert( ~isempty(fname), 'Missing path to config file (env:OSLCONF).' );

			assert( isstruct(s) && isscalar(s), 'First input to osl_conf.write() must be a struct.' );
			assert( all(structfun( @ischar, s )), 'Config values must be strings.' );
			
			% overwrite config file
			txt = structkvfun( @(k,v) sprintf('%s = %s',k,v), s, false );
			filewrite( fname, txt );

		end

	end

end