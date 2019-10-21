function osl_upgrade
    % This function calls osl_core/upgrade.sh to update a non-Git release
    % with the latest code from GitHub
    %
    % WARNING - this function MAY break your OSL release. Use with caution
    %
    system(sprintf('cd %s && ./upgrade.sh',fullfile(osldir,'osl-core')))