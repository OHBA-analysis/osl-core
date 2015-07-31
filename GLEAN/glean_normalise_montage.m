function D = glean_normalise_montage(S)
% Obsolete code for performing normalisation

D = spm_eeg_load(S.D);

samples2use = find(~all(badsamples(D,':',':',':')));

if strcmp(D.transformtype,'TF')
    stdev = std(D(:,:,samples2use,1),[],3);
else
    stdev = std(D(:,samples2use,1),[],2);
end
switch(S.normalisation)

    case 'voxel'
        tra  = diag(1./mean(stdev,2));
        name = 'voxelwise normalisation';
        
    case 'global' % using mean(std(x)) instead of std(x(:))
        tra = eye(D.nchannels)./mean(stdev,2);
        name = 'global normalisation';
        
    case 'none'
        tra = eye(D.nchannels);
        name = 'no normalisation';
        
end
        
% Write normalisation as a montage:
montage = [];
montage.labelorg = D.chanlabels';
montage.labelnew = D.chanlabels';

montage.name = name;
montage.tra = tra;
montage.chantypenew = D.chantype';
montage.chantypeorg = D.chantype';


S2 = [];
S2.D            = fullfile(D.path,D.fname);
S2.montage      = montage;
S2.keepsensors  = false;
S2.mode         = 'write';

D = spm_eeg_montage(S2);
D.save;

end