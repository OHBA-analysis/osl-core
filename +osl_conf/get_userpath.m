function upath = get_userpath()
%
% Locate the directory to read/write OSL configurations to/from
%
% This is (in order of priority) either:
% 1) set by the environment variable OSLUSERDIR
% 2) the directory of the osl.configuration file
%
% Create a default config file in OSLDIR otherwise.
%
% JH

    f = getenv('OSLUSERDIR'); % manually set
    
    if ~isempty(f)
        assert( osl_util.isdir(f), 'OSL user-path not found (set via $OSLCONF): %s', f );
        return;
    end


end
