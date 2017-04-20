function osl_upgrade
	% This function calls osl_core/upgrade.sh to update a non-Git release
	% with the latest code from GitHub
	system(sprintf('cd %s && ./upgrade.sh',fullfile(osldir,'osl-core')))