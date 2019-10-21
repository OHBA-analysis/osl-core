function y = isfile( fname, varargin )
%
% ans = osl_util.isfile( fileName )
%
% Complement Matlab's function isdir().
% Returns logical for exist(fname,'file') = [2,3,4,6]
%
% See also: osl_util.isdir
%
% JH

    % NOTE:
    % Keep varargin for backwards compat when converting from 
    %   exist( filename, 'file' ) == 2    =>    osl_util.isfile(filename, 'file')

    y = [false, false, true, true, true, false, true, false, false];
    y = y(1 + exist(fname,'file'));

end