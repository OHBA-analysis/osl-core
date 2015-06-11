function [data,freq_ind] = bandpass(data,freqrange,fsample,do_plot,freq_ind)
% data = bandpass(data,fbands,fsample,do_plot)
% note freq_ind can also be passed in if already computed from a previous
% call.
% data - Nsignals x Ntpts

if any(size(data)==1)
    data = data(:)';
end

if(nargin < 4)
    do_plot = 0;
end

t = (1:size(data,2))/fsample;

if do_plot
    h1 = figure;
    h2 = figure; 
    plot(t,data,'r')
end

data_fft = fft(data,[],2);

fHz = linspace(0,fsample,size(data,2));

if nargin > 4 && ~isempty(freq_ind)
    bp_filter = zeros(size(fHz));
    bp_filter(freq_ind) = 1;
else
    bp_filter = fHz > freqrange(1) & fHz <= freqrange(2);
    freq_ind = find(bp_filter);
    bp_filter = bp_filter | fliplr(bp_filter);
end

data_fft(:,~bp_filter) = 0;

data = real(ifft(data_fft,[],2));

if do_plot
    figure(h1)
    plot(fHz(1:ceil(size(data_fft,2)/2)),abs(data_fft(:,1:ceil(size(data_fft,2)/2))))
    figure(h2), hold on
    plot(t,data,'b')
end

end
