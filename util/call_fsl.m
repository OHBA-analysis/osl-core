function [status,output] = call_fsl(cmd)
	% Deprecated, left in for compatibility
	% Note how the order of the outputs is swapped to maintain compatibility
	fprintf(2,'Warning - this function is deprecated, preferred usage is to call runcmd() directly\n');
	[output,status] = runcmd(cmd);
