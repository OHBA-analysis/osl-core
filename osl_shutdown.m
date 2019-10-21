function osl_shutdown(force)
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

    if nargin < 1, force=false; end

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
    if force
        p = strsplit(path,pathsep);
        p = p(cellfun( @(x) osl_util.contains(x,osl_root), p ));
        cellfun( @(x) rmpath(x), p );
    end

    % Done - clear OSLDIR
    setenv('OSLDIR');
    fprintf(1,'[OSL] Cleaned up the path. Bye now!\n');

end

