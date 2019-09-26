function bf_save(BF, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_save.m 104 2014-01-15 17:45:40Z litvak.vladimir@gmail.com $
%
% This is a patch that allows us saving in -v7, otherwise incredibly bad performance 
% System: MBP '15, Matlab R2016a, OSx 10.11.6 
% RB 2017

    if nargin < 2, overwrite=true; end

    fpath = fullfile(pwd,'BF.mat');
    if osl_util.isfile(fpath) && ~overwrite
        try 
            save(fpath, '-struct', 'BF', '-append');
        catch
            fprintf(2,'Could not append to file: %s\n',fpath);
            fprintf(2,'This could be because the storage format is v7, and appending would make the file larger than 2GB.\n');
            fprintf(2,'In that case, either manually load and re-save with format v7.3, or set overwrite to true.\n');
        end
    else
        % appending to -v7.3 slows down performance incredibly, so modified here
        % todo: is that sufficient storage capacity for us?

        % JH: changed here to adapt to larger file-sizes
        bs = osl_util.bytesize(BF);
        if bs.GB >= 2
            save(fpath, '-struct', 'BF', '-v7.3');
        else
            save(fpath, '-struct', 'BF', '-v7');
        end 
    end
    
end
