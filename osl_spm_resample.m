function resampled_signal = osl_spm_resample(signal,ds_factor)
% this is to downsample the time domain signal.  It's the spm
% function 'spm_resample'

if(size(signal,2)==1)
    error('can not downsample singleton dimension');
end;

N0          = size(signal,2);
N           = floor(N0*ds_factor);
ds_factor   = N/N0;
Y           = fftshift(fft(signal,[],2),2);
sy          = size(Y,2);
middle      = floor(sy./2)+1;

if ds_factor>1 % upsample
    N2 = floor((N-N0)./2);
    if N0/2 == floor(N0/2)
        Y(:,1) = []; % throw away non symmetric DFT coef
    end
    Y  = [zeros(size(Y,1),N2),Y,zeros(size(Y,1),N2)];
else % downsample
    N2 = floor(N./2);
    Y  = Y(:,middle-N2:middle+N2,:,:);
end
resampled_signal = ds_factor*ifft(ifftshift(Y,2),[],2);
