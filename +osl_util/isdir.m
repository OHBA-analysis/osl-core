function y = isdir( dname, varargin )
%
% ans = osl_util.isdir( folderName )
%
% Alternative to Matlab's isdir function, using exist.
%
% See also: osl_util.isfile
%
% JH

    y = exist(dname,'dir') == 7;

end