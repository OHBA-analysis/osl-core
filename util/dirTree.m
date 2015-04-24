function filelist = dirTree(dirname)
% filelist = dirTree(dirname)
% Finds all files in directory DIRNAME by looking in each subfolder and
% returns a list of files FILELIST
%
% Adam Baker 2015

filelist = [];
directory = dir(dirname);
for i = 1:length(directory)
    if(~strcmp(directory(i).name,'.') && ~strcmp(directory(i).name,'..'))
        if(directory(i).isdir)
            filelist = [filelist; dirTree(fullfile(dirname,directory(i).name))];
        else
            filelist = [filelist; {fullfile(dirname,directory(i).name)}];
        end
    end
end
end