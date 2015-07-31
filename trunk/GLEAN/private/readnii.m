function mat = readnii(nii_fname,mask_fname)

mat = vols2matrix(read_avw(nii_fname), read_avw(mask_fname));

end