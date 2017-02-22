function mat = readnii(nii_fname,mask_fname)

if nargin > 1 
    mat = vols2matrix(read_avw(nii_fname), read_avw(mask_fname));
else
    mat = read_avw(nii_fname);
end

end