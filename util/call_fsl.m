function [status,output] = call_fsl(cmd)
	% OSL handles FSL setup, so this function is used to override the
	% version of call_fsl.m that comes with FSL to avoid issues related
	% to environment variables being cleared in the Matlab session
	%
	% Note how the order of the outputs is swapped compared to runcmd() to maintain compatibility
	% with other FSL matlab tools (like nii.load)
	[output,status] = runcmd(cmd);
