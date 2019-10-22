function s = osldir(varargin)
% Return the OSL root directory stored in environment variable 'OSLDIR'

    s = getenv('OSLDIR');
    assert( ~isempty(s), 'Environment variable OSLDIR is not set.' );
    s = fullfile( s, varargin{:} );
    
end