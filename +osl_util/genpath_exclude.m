function pathstr = genpath_exclude(pathstr,excludes)
    % Take in list of strings to exclude from path

    if ischar(excludes)
        excludes = {excludes};
    end

    paths = genpath(pathstr);
    paths = strsplit(paths,pathsep);
    retain = true(size(paths));

    for j = 1:length(excludes)
        retain = retain & cellfun( @(x) ~osl_util.contains(x,excludes{j}), paths );
    end

    paths = paths(retain);
    pathstr = strjoin(paths,pathsep);

end