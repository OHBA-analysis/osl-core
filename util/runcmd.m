function [output,return_value] = runcmd(cmd,varargin)
	% Wrapper to capture output of system calls. Supports sprintf token replacement
	% 
	% USAGE
	% - runcmd('ls')
	% - output = runcmd('ls')
	% - output = runcmd('ls %s',getenv('OSLDIR'))


	if nargin > 1
		cmd = sprintf(cmd,varargin{:});
	end

	[return_value, output] = system(cmd);

	if(return_value ~= 0)
	    throw(MException('runcmd:error',sprintf('runcmd call:\n%s\nReturn value: %d\nProduced error:\n%s\n',cmd,return_value,output)));
	end
