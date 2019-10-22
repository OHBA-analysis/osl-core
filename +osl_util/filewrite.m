function count = filewrite( fpath, txt, fmode )
%
% count = osl_util.filewrite( fpath, txt, mode='w+' )
%
% Complement Matlab's function fileread.
% This function writes or appends text `txt` to file `fpath`.
%
% fpath can be an absolute or relative path, or just a filename (will be saved in current directory).
% txt should be a string, or a cell of strings, in which case cells are concatenated as individual lines.
%
% mode determines what happens if the file already exists:
%   w      Overwrite contents of already existing file (error if not found).
%   w+     DEFAULT: Overwrite contents of file, create it if needed. 
%   a      Append to contents of existing file (error if not found).
%   a+     Append to file if it exists, otherwise create it.
%
% Output count corresponds to the output of fwrite (number of chars written).
%
% See also: fopen, fwrite
%
% JH

    % process fmode
    if nargin < 3, fmode='w+'; end
    assert( ismember(fmode, {'w','a','w+','a+'}), 'Invalid mode: %s', fmode );

    % process text
    if iscellstr(txt)
        txt = strjoin(txt,newline);
    end
    assert( ischar(txt), 'Input text should be a string or a cell of strings.' );

    % make sure parent folder exists
    folder = fileparts(fpath);
    assert( isempty(folder) || isdir(folder), 'Folder not found: %s', folder );

    % write to file
    fh = fopen(fpath, fmode);
    count = fwrite(fh, txt);
    fclose(fh);

end
