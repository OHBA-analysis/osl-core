function D = glean_normalise(S)
% Performs normalisation of MEEG data

D = spm_eeg_load(S.D);

% TODO: Gotta sort this out for trialwise data
trl = 1;

nFreq = D.nfrequencies;
if isempty(nFreq)
    nFreq = 1;
end

% Maybe block this:
if isequal(D.transformtype,'TF')
    means = mean(D(:,:,:,:),3);
    stdev = std(D(:,:,:,:),[],3);
else
    means = mean(D(:,:,:),2);
    stdev = std(D(:,:,:),[],2);
end
            
            
for f = 1:nFreq
    
    if isequal(D.transformtype,'TF')
        
        D(:,f,:,trl) = D(:,f,:,trl) - repmat(means(:,f),1,1,D.nsamples);
        switch(S.normalisation)
            case 'voxel'
                D(:,f,:,trl) = D(:,f,:,trl) ./ repmat(stdev(:,f),1,1,D.nsamples);
            case 'global'
                D(:,f,:,trl) = D(:,f,:,trl) ./ repmat(mean(stdev(:,f)),D.nchannels,1,D.nsamples);
        end
        
    else
        
        D(:,:,trl) = D(:,:,trl) - repmat(means,1,D.nsamples);
        switch(S.normalisation)
            case 'voxel'
                D(:,f,:,trl) = D(:,:,trl) ./ repmat(stdev,1,D.nsamples);
            case 'global'
                D(:,f,:,trl) = D(:,:,trl) ./ repmat(mean(stdev),D.nchannels,D.nsamples);
        end
        
    end
    
end

D.save;


end