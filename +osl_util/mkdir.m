function mkdir( dname )
%
% osl_util.mkdir( folderName )
% 
% Call Matlab's mkdir only if folder does not exist (avoiding errors if it does).
%
% JH

    if ~osl_util.isdir(dname)
        mkdir(dname); % will error on failure
    end

end