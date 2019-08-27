function y = isfilefname, varargin)
%
%  res = osl_util.isfile( fileName )
%
% Complement Matlab's function isdir().
% Returns logical for exist(fname,'file') = [2,3,4,6]
%
% JH

    % NOTE:
    % Keep varargin to prevent silly error when converting from 
    %   exist( filename, 'file' ) == 2    =>    osl_util.isfile(filename, 'file')

    y = [false, false, true, true, true, false, true, false, false];
    y = y(1 + exist(fname,'file'));

end
