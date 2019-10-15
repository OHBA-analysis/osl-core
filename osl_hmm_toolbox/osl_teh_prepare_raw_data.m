function [data, lags] = osl_teh_prepare_raw_data(data,normalisation,logtrans,embed)

lags=[];

% log transform
if logtrans
    data = log10(data);
end

lags=[];

if exist('embed','var') && embed.do
    
    disp('Time embedding data');
    
    if isfield(embed,'centre_freq')
        span=1/embed.centre_freq; %secs    
        num_embeddings=round(span/embed.tres); 
    elseif isfield(embed,'num_embeddings')
        num_embeddings=embed.num_embeddings; 
    end
    
    %lags=round(linspace(-num_embeddings/2,num_embeddings/2,num_embeddings));
    lags=floor(-num_embeddings/2:num_embeddings/2);
    %lags=round(0:num_embeddings-1);

    %disp(lags);
    
    [dataout,valid]=embedx(data',lags);
    dataout=dataout';
    data=randn(size(dataout,1),size(data,2))*std(squash(data(:,1:500)));
    data(:,valid)=dataout;          
       
end

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
    end
end

end