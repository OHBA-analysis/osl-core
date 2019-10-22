function y = generate_ar_surrogates(data,ar_order)
    % Generate surrogate data by fitting an AR model using the Yule-Walker method
    % to each channel independently, and then generating samples by randomly driving
    % each channel. As channels are independent, there should be no genuine correlations
    % between channels in the surrogate

    if nargin < 2 || isempty(ar_order)
        error('Must specify the order of the AR model');
    end
    
    y = nan(size(data));
    means = mean(data);
    data = bsxfun(@minus,data,means);

    for j = 1:size(data,2)
        R=xcorr(data(:,j),ar_order,'biased');
        R = R(ar_order+1:end);
        cf = inv(toeplitz(R(1:end-1)))*R(2:end); % AR coefficients
        stdev = sqrt(R(1) - R(2:end)'*cf);

        V = stdev*randn(size(data,1),1);
        y(:,j) = filter(1,[1;-cf(:)],V);

    end

    y = bsxfun(@plus,y,means);
