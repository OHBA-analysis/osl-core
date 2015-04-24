function sortedList = osl_filelist(pathstr,filestr)
%OSL_FILELIST Gets a list of files from a directory matching a pattern
%
% LIST = OSL_FILELIST(DIRNAME, PATTERN) matches files with string PATTERN
%   in directory DIRNAME, sorting in a sensible order. 
%
%
% Adam Baker 2014


if iscellstr(pathstr)
    pathstr = pathstr{:};
elseif ~ischar(pathstr)
    error('pathstr must be a string')
end

if ~strcmp(pathstr(end),filesep)
    pathstr = strcat(pathstr,filesep);
end
tmp        = dir(fullfile(pathstr,filestr));
filelist   = {tmp(:).name};
filelist   = cellfun(@(x) strcat(pathstr,x),filelist,'uniformoutput',0);
sortedList = sort_nat(filelist);
sortedList = sortedList(:);
    
end