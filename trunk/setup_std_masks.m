function setup_std_masks


gridsteps=[2:15];

for i=1:length(gridsteps),
    i
    gridstep=gridsteps(i);
    stdmasksdir='/home/mwoolrich/vols_data/std_masks';

    runcmd(['fslcreatehd ' num2str(round([91 109 91]*2/gridstep)) ' 1 ' num2str(gridstep) ' ' num2str(gridstep) ' ' num2str(gridstep) ' 1 0 0 0 4 ' stdmasksdir '/MNI152_T1_brain_tmp.nii.gz']); 

    runcmd(['flirt -in ' getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain -applyxfm -init ' getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' stdmasksdir '/MNI152_T1_' num2str(gridstep) 'mm_brain -paddingsize 0.0 -interp trilinear -ref ' stdmasksdir '/MNI152_T1_brain_tmp']);

    runcmd(['flirt -in ' getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask -applyxfm -init ' getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' stdmasksdir '/MNI152_T1_' num2str(gridstep) 'mm_brain_mask -paddingsize 0.0 -interp nearestneighbour -ref ' stdmasksdir '/MNI152_T1_brain_tmp']);
    
    runcmd(['fslmaths ' stdmasksdir '/MNI152_T1_' num2str(gridstep) 'mm_brain -mas ' stdmasksdir '/MNI152_T1_' num2str(gridstep) 'mm_brain_mask ' stdmasksdir '/MNI152_T1_' num2str(gridstep) 'mm_brain']);

end;
