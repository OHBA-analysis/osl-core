function [data, lags] = osl_teh_prepare_data(D,normalisation,logtrans,freq_ind,embed)

% [data, lags] = osl_teh_prepare_data(D,normalisation,logtrans,freq_ind,embed)
%
% prepares data (including time embedding) for calling HMM-MAR

% reshape trialwise data
if strcmp(D.transformtype,'TF')
    data = D(:,freq_ind,:,:);
else
    data = D(:,:,:);
end
data = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);

% select only good data
goodsamples = good_samples(D);
goodsamples = reshape(goodsamples,1,D.nsamples*D.ntrials);
data = data(:,goodsamples);

[data, lags] = osl_teh_prepare_raw_data(data,normalisation,logtrans,embed);

