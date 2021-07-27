function p = path_backup()
%
% Used at startup to backup the Matlab path.
% When calling osl_shutdown, the backup is used to restore the original path.
% 
% See also: osl_startup, osl_util.path_restore
%
% JH

    p = path(); 

    % always backup next to current OSLCONF
    confdir = fileparts(getenv('OSLCONF'));
    fprintf(1,'[OSL] Backing-up matlab path to: %s\n',confdir);

    fname = fullfile( confdir, '.path-backup.tmp' );
    setenv('OSL_PATH_BACKUP', fname);
    fh = fopen(fname,'w+');
    fwrite(fh,p);
    fclose(fh);

    p = strsplit(p,pathsep);

end
