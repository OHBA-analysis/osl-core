function [ret, w]=runcmd(cmd,varargin)
	% Wrapper for system calls
	% Added functionality - cmd is passed through sprintf so string substitutions can be
	% applied in runcmd() directly
	% EXAMPLE USAGE
	% runcmd('rm %s','testdir')

	if nargin > 1
		cmd = sprintf(cmd,varargin{:});
	end

	[ret, w] = system(cmd);

	if(ret ~= 0),
	    throw(MException('runcmd:error',sprintf('runcmd call:\n%s\nReturn value: %d\nProduced error:\n%s\n',cmd,ret,w)));
	end;
