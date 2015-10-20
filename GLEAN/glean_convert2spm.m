function D = glean_convert2spm(filename,data,fs)
% Converts data to SPM12 format, saving the new .mat and .dat files.
%
% D = GLEAN_CONVERT2SPM(filename,data,fs)
%
% REQUIRED INPUTS:
%   filename  - Name of the SPM12 MEEG object to create
%   data      - [channels x samples] data to convert
%   fs        - sampling rate of the data in Hz
% 
% OUTPUTS:
%   D         - Newly created SPM12 MEEG object 
%
% Adam Baker 2015


[pathstr,filestr] = fileparts(filename);
if isempty(pathstr)
    pathstr = pwd;
end
filename = fullfile(pathstr,filestr);

blank = fullfile(fileparts(mfilename('fullpath')),'data','blank_meeg');
D = spm_eeg_load(blank);
D = clone(D,filename,[size(data,1),size(data,2),size(data,3)]);
D(:,:,:) = data;
D = fsample(D,fs);
D.save;

end