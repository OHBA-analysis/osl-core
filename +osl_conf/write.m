function write(s,fname)
%
% osl_conf.write( s, fname = $OSLCONF )
%
% Write config structure to file.
%
% JH

    if nargin < 2, fname=getenv('OSLCONF'); end
    assert( ~isempty(fname), 'Missing path to config file.' );

    if ~isempty(getCurrentTask)
        fprintf(2,'Running on parallel worker, not writing to "%s" to avoid corruption.\n',fname);
        return
    end
    
    assert( isstruct(s) && isscalar(s), 'Config must be a struct.' );
    assert( all(structfun( @ischar, s )), 'Config values must be strings.' );
    
    txt = osl_util.structkvfun( @(k,v) sprintf('%s = %s',k,v), s, false );
    osl_util.filewrite( fname, txt );
    
end