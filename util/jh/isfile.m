function y = isfile(fname)
%
%  res = isfile( fileName )
%
% Complement Matlab's function isdir().
% Returns logical for exist(fname,'file') = [2,3,4,6]
%
% JH

    y = [false, false, true, true, true, false, true, false, false];
    y = y(1 + exist(fname,'file'));

end