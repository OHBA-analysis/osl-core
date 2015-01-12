function filelist = osl_filelist(pathstr,filestr)
% Gets a list of files from a directory specified by pathstr that match the
% in filestr
%
% Usage:
% filelist = osl_filelist(dirname,filestr)
%
%
% Adam Baker 2014
    filelist = dir(fullfile(pathstr,filestr));
    filelist = {filelist.name};
    filelist = fullfile(pathstr,filelist);
    filelist = sort_nat(filelist)';
    
end