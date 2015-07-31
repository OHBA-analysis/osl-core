function data_out = dataio(filename,data_in)
% data_out = DATAIO(filename,data_in)
%
% Write data to or read data from binary file. For large arrays this is 
% faster than saving as a .mat file
%
% USAGE: 
% -------------------------------------------
% DATAIO(filename,data) WRITES the N dimensional array DATA
% to the file specified by FILENAME
%
% data = DATAIO(filename) READS data from the file specified by FILENAME
% into an N dimensional array DATA
% -------------------------------------------
% AB 2012



% global settings:
maxdims    = 5; % maximum dimensionality
maxdimsize = 8; % maximum dimension size allowed 

% do some checks
if nargin == 0
  error('function requires filename for read mode and filename and data for write mode')
elseif nargin == 1
  mode = 'r';
elseif nargin == 2 && nargout == 0
  mode = 'w';
end

[~, ~, ext] = fileparts(filename);
if strcmp(ext,'')
  filename = [filename '.bin'];
elseif ~strcmp(ext,'.bin')
  error('filename must be a binary file')
end

switch mode
  case 'r'
    data_out = readdat(filename);
  case 'w'
    savedat(filename,data_in);
end

   

function savedat(filename,data)

% get dimensionality of data:
datasize_new = size(data);
if exist(filename,'file')
    datasize_org = getdatasize(filename);
else
    datasize_org = 0;
end

if all(datasize_new(1:end-1) == datasize_org(1:end-1))
    datasize = [datasize_new(1:end-1), datasize_org(end) + datasize_new(end)];
else
    error('Incompatible data sizes');
end

% check if filename exists
[pathstr, ~, ~] = fileparts(filename);
if ~isdir(pathstr)
  mkdir(pathstr);
end


% create header:
hdr = sprintf(sprintf('%%0%02ii',maxdims*maxdimsize),0);
for dim = 1:length(datasize)
  str = sprintf(sprintf('%%%02ii',maxdimsize),datasize(dim));
  hdr((1:maxdimsize) + (maxdimsize*(dim-1))) = str;
end


% write header to file:
if ~exist(filename,'file')
    fid = fopen(filename,'w');
    fwrite(fid,hdr,'char');
    fclose(fid);
else
    fid = fopen(filename,'r+');
    frewind(fid);
    fwrite(fid,hdr,'char');
    fclose(fid);
    
end

fid = fopen(filename,'a');
fwrite(fid,data,'double');
fclose(fid);

end



function data = readdat(filename)

% read in header and data
fid = fopen(filename,'r');
hdr = fread(fid,[1,40],'*char');
datasize = str2num(reshape(hdr,[8,5])')';
datasize(datasize==0) = [];
data = fread(fid,datasize,'*double');
%data = reshape(data,datasize(end:-1:1));
fclose(fid);

end


function datasize = getdatasize(filename)

% read in header and data
fid = fopen(filename,'r');
hdr = fread(fid,[1,40],'*char');
datasize = str2num(reshape(hdr,[8,5])')';
datasize(datasize==0) = [];
fclose(fid);

end


end