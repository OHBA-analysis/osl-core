function writenii(map,fname,mask_fname)
% Write data to nifti file, based on an input mask

dirname = fileparts(fname);
if ~isdir(dirname)
    mkdir(dirname);
end

[mask,~,scales] = read_avw(mask_fname);
save_avw(matrix2vols(map,mask),fname,'f',scales);

gunzip([fname '.gz']);
system(['rm ' fname '.gz']);

system(['fslcpgeom ' mask_fname ' ' fname ' -d']);

end