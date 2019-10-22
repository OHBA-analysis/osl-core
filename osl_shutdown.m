function osl_shutdown( clean_backup_path )
% 
% osl_shutdown( clean_backup_path=false )
% 
% Restore Matlab path to its original value (before osl_startup was called).
%
% The folder osl-core itself is often manually added to the path prior to 
% calling osl_startup (in order to make osl_startup callable). Hence, the 
% restored path will often still contain osl-core. Set clean_backup_path to
% true in order to remove OSL folders from the backed up path as well.
%
% RA 2017, JH 2019

    if nargin < 1, clean_backup_path=false; end

    % OSL is initialized if OSLDIR is set. If OSLDIR is not set, then do nothing
    if osl_isactive()
        osl_root = getenv('OSLDIR');
    else
        return
    end

    % Restore path backup
    % JH: use separate file for path backup
    osl_conf.path_restore();

    % It's concievable that the backup path still contains OSL directories; for example
    % it is common to add the osl-core folder manually in order to call osl_startup.
    if clean_backup_path
        p = strsplit(path,pathsep);
        p = p(cellfun( @(x) ~isempty(strfind(x,osl_root)), p )); %#ok
        cellfun( @(x) rmpath(x), p );
    end

    % Done - clear OSLDIR
    setenv('OSLDIR');
    fprintf(1,'[OSL] Cleaned up the path. Bye now!\n');

end

