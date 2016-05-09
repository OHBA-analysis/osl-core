function [data, lags] = osl_teh_prepare_data(D,normalisation,logtrans,freq_ind,embed)

% [data, lags] = osl_teh_prepare_data(D,normalisation,logtrans,freq_ind,embed)
%
% prepares data (including time embedding) for calling HMM-MAR

lags=[];

% reshape trialwise data
if strcmp(D.transformtype,'TF')
    data = D(:,freq_ind,:,:);
else
    data = D(:,:,:);
end
data = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);

% select only good data
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
data = data(:,good_samples);

% log transform
if logtrans
    data = log10(data);
end

if exist('embed','var') && embed.do,    
    
    disp('Time embedding data');
    span=1/embed.centre_freq; %secs
    num_embeddings=round(span/embed.tres);        
    %lags=round(linspace(-num_embeddings/2,num_embeddings/2,num_embeddings));
    lags=round(-num_embeddings/2:num_embeddings/2);
    %lags=round(0:num_embeddings-1);

    %disp(lags);
    
    [dataout,valid]=embedx(data',lags);
    dataout=dataout';
    data=randn(size(dataout,1),size(data,2))*std(squash(data(:,1:500)));
    data(:,valid)=dataout;          
       
end;

% normalisation
switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
end

if embed.rectify
    data=abs(data);
    
    switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
    end;
end

end
