
function filestr = prefix(filestr, prefixstr)
% PREFIX prefixes a character onto a set of filenames.
% FILESTR = PREFIX(FILESTR, PREFIXSTR) adds prefix PREFIXSTR to filename(s)
% to filenames in FILESTR. FILESTR can be a single string, or a cell
% array of strings. PREFIXSTR is a character string.

assert(ischar(prefixstr),               ...
    [mfilename ':IncorrectPrefixType'], ...
    'Prefix should be a string. \n');

% apply prefix, catering for character or cell arrays
if iscell(filestr)
    for s = 1:numel(filestr)
        [p,f,e] = fileparts(filestr{s});
        filestr{s} = fullfile(p,[prefixstr,f,e]);
    end
elseif ischar(filestr)
    [p,f,e] = fileparts(filestr);
    filestr = fullfile(p,[prefixstr,f,e]);
else
    error([mfilename ':InconsistentInput'], ...
        'Expecting a string or cell array of strings as input. \n');
end


end