function osl_shutdown
% Remove OSL from path
%
% Note that this will not remove FSL or workbench from the system path
% so commands like "system('fslval')" should still run
% In that respect, there will still be some traces of OSL. It's a bit more 
% complicated to remove these because they may not have been added by
% OSL, so best not to worry about them for now. The FSL Matlab helpers
% may also be left intact, if they were present before OSL was started
%
% Romesh Abeysuriya 2017
% JH 2018

	% OSL is initialized if OSLDIR is set. If OSLDIR is not set, then do nothing
	if isempty(getenv('OSLDIR'))
		return
	else
		osl_root = getenv('OSLDIR');
	end

	% Restore path backup
	% JH: use separate file for path backup
	restore_path();

	% It's concievable that the backup path still contains OSL directories; for example
	% it is common to add the osl-core folder manually in order to call osl_startup.
	p = strsplit(path,pathsep);
	p = p(cellfun( @(x) ~isempty(strfind(x,osl_root)), p ));
	cellfun( @(x) rmpath(x), p );

	% Done - clear OSLDIR
	setenv('OSLDIR');
	fprintf(1,'[OSL] Cleaned up the path. Bye now!\n');

end

function restore_path()

	fname = getenv('OSL_PATH_BACKUP');
	if ~isempty(fname)
		if exist(fname,'file') ~= 2
			fprintf(2,'Unable to restore path due to missing file: %s\n',fname);
		else
			p = fileread(fname);
			if ~isempty(p)
				path(p);
            end
            %delete(fname);
		end
	end
	setenv('OSL_PATH_BACKUP');

end