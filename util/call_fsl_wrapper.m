function result = call_fsl_wrapper(fslCommand, quiet)
	% Deprecated, left in for compatibility
	% Adds option to echo command to terminal

	if nargin < 2 || ~exist('quiet', 'var') || ~quiet,
	    fprintf(fslCommand);
	    fprintf('\n');
	end

	[status, result] = call_fsl(fslCommand);
