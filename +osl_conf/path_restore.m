function path_restore()
%
% Used at shutdown to restore original Matlab path.
% 
% See also: osl_shutdown, osl_util.path_backup
%
% JH

    fname = getenv('OSL_PATH_BACKUP');
    if isempty(fname)
        warning( '[bug] Path backup is not set.' );
    else
        if ~osl_util.isfile(fname)
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