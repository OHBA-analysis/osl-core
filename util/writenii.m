function writenii(map,fname,mask_fname)
% Write data to nifti file, based on an input mask

[dirname,~,~] = fileparts(fname);
if ~isdir(dirname)
    mkdir(dirname);
end

if nargin > 2
    % Save with geometry
    [mask,~,scales] = read_avw(mask_fname);
    save_avw(matrix2vols(map,mask),fname,'f',scales);
    system(['fslcpgeom ' mask_fname ' ' fname ' -d']);
else
    % Save as size of map
    save_avw(map,fname,'f',[1 1 1 1]);
end
    

end
