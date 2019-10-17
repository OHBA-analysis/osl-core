function v = osl_version
    % Read and return the version number for the OSL distribution from the version.txt file
    % If this file doesn't exist, then return the commit hash for osl-core
    % If git is not functional, return 'unknown'
    
    version_file = fullfile(osldir,'version.txt');
    if osl_util.isfile(version_file)
        fid = fopen(version_file,'r');
        l = fgetl(fid);
        fclose(fid);
        v = strrep(l,' = ','');
    else
        try
            v = runcmd('git -C %s rev-parse HEAD ',fullfile(osldir,'osl-core'));
            v = v(1:7);
        catch
            v = 'unknown';
        end
    end
