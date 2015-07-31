function newNames = prefix(oldNames, prefix)
%PREFIX  prefixes a character onto a set of filenames. 
% NEW = PREFIX(OLD, PRE) adds prefix PRE to filename OLD to make
%   NEW filename. OLD can be a single string, or a cell array of strings.
%   PRE is a character string. OLD can contain directories, these are taken
%   care of.

% input checking
assert(ischar(prefix),                     ...
       [mfilename ':IncorrectPrefixType'], ...
       'Prefix should be a string. \n');
   
% apply prefix, catering for character or cell arrays
if iscell(oldNames), 
    newNames = cellfun(@(name) add_prefix(name, prefix), ...
                       oldNames,                         ...
                       'UniformOutput', false);
    
elseif ischar(oldNames),
    newNames = add_prefix(oldNames, prefix);
    
else 
    error([mfilename ':InconsistentInput'], ...
          'Expecting a string or cell array of strings as input. \n');
end%if
end%prefix

function new = add_prefix(old, pre)
%ADD_PREFIX

[dirName, fName, ext] = fileparts(old);

new = fullfile(dirName, sprintf('%s%s%s', pre, fName, ext));
end%add_prefix
% [EOF]