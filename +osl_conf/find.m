function f = find()
%
% Locate the file osl.conf.
%
% The order of priority is:
%  - Manually set via environment variable OSLCONF;
%  - In userpath if $OSLUSERMODE is 'user'
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
        % Config specified and not found - error
        assert( osl_util.isfile(f), 'Configuration file not found (set via $OSLCONF): %s', f );
        return;
    elseif isempty(f) && strcmp(getenv('OSLUSERMODE'), 'user')
        % Config not specified and in user mode - set fname to userpath
        f = fullfile(userpath, 'osl.conf');
    elseif osl_util.isfile(p)
        % Config not specified use $OSLDIR
        f = p;
    elseif osl_util.isfile(c)
        % Config not found in $OSL - use osl-core directory
        f = c;
    end

    if ~osl_util.isfile(f)
        % no config found in specified location - create a default
        s = osl_conf.default();

        warning( 'Configuration file osl.conf not found; creating default "%s".', f );
        osl_conf.write(s,f);
    end

end
