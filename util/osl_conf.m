classdef (Abstract) osl_conf

	methods(Static)

		function s = read()
	
			fname = getenv('OSLCONF');
			assert( ~isempty(fname), 'Missing path to config file (env:OSLCONF).' );

			if ~osl_util.isfile(fname)
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

			fname = getenv('OSLCONF');
			assert( ~isempty(fname), 'Missing path to config file (env:OSLCONF).' );
            
            if ~isempty(getCurrentTask)
				fprintf(2,'Running on parallel worker, not writing to "%s" to avoid corruption.\n',fname);
				return
			end

			assert( isstruct(s) && isscalar(s), 'First input to osl_conf.write() must be a struct.' );
			assert( all(structfun( @ischar, s )), 'Config values must be strings.' );
			
			% overwrite config file
			txt = osl_util.structkvfun( @(k,v) sprintf('%s = %s',k,v), s, false );
			filewrite( fname, txt );

		end

	end

end