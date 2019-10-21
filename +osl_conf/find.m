function f = find()
%
% Locate the file osl.conf. 
%
% The order of priority is:
%  - Manually set via environment variable OSLCONF;
%  - In $OSLDIR;
%  - In osl-core directory.
%
% Create a default config file in OSLDIR otherwise.
%
% JH

    p = osldir('osl.conf');
    c = fileparts(mfilename('fullpath'));
    c = fullfile( c, '..', 'osl.conf' );
    f = getenv('OSLCONF'); % manually set
    
    if ~isempty(f)
        assert( osl_util.isfile(f), 'Configuration file not found (set via $OSLCONF): %s', f );
        return;
    end

    if osl_util.isfile(p)
        f = p;
    elseif osl_util.isfile(c)
        f = c;
    else
        s = osl_conf.default();
        f = p;
        
        warning( 'Configuration file osl.conf not found; creating default "%s".', f );
        osl_conf.write(s,f);
    end

end