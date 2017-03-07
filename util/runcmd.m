function [output,return_value] = runcmd(cmd,varargin)
	% Wrapper for system calls
	% Added functionality - cmd is passed through sprintf so string substitutions can be
	% applied in runcmd() directly
	% EXAMPLE USAGE
	% runcmd('rm %s','testdir')

	if nargin > 1
		cmd = sprintf(cmd,varargin{:});
	end

	[return_value, output] = system(cmd);

	if(return_value ~= 0)
	    throw(MException('runcmd:error',sprintf('runcmd call:\n%s\nReturn value: %d\nProduced error:\n%s\n',cmd,return_value,output)));
	end
