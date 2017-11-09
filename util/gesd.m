function [idx,x2] = gesd(x,alpha,n_out,outlier_side)
    % Detect outliers using Generalized ESD test
    %
    % INPUTS
    % - x : Data set containing outliers
    % - alpha : Significance level to detect at (default = 0.05)
    % - n_out : Maximum number of outliers to detect (default = 10% of data set)
    % - outlier_side : Specify sidedness of the test
    %       - outlier_side = -1 -> outliers are all smaller
    %       - outlier_side = 0  -> outliers could be small/negative or large/positive (default)
    %       - outlier_side = 1  -> outliers are all larger
    % 
    % OUTPUTS
    % - idx : Logical array with TRUE wherever a sample is an outlier
    % - x2 : Input array with outliers removed
    %
    % For details about the method, see
    % B. Rosner (1983). Percentage Points for a Generalized ESD Many-Outlier Procedure. Technometrics 25(2), pp. 165-172. 
    % http://www.jstor.org/stable/1268549?seq=1
    %
    % Romesh Abeysuriya 2013-2017
    
    if nargin < 4 || isempty(n_out)
        outlier_side = 0; % By default the outliers can be in both directions
    end

    if nargin < 3 || isempty(n_out)
        n_out = ceil(length(x)*0.1); % Assume no more than 10% of the data is outliers by default
    end

    if nargin < 2 || isempty(alpha)
        alpha = 0.05;
    end

    if outlier_side == 0;
        alpha = alpha/2;
    end

    if any(~isfinite(x))
        % Need to find outliers only in non-finite x
        y = find(isfinite(x)); % These are the indexes of x that are finite
        [idx1,x2] = gesd(x(isfinite(x)),alpha,n_out,outlier_side);
        % idx1 has the indexes of y which were marked as outliers
        % the value of y contains the corresponding indexes of x that are outliers
        idx = false(1,length(x));
        idx(y(idx1)) = true;
        return
    end        

    n = length(x);
    temp = x;
    R = zeros(1,n_out);
    rm_idx = R;
    lambda = R;
    
    for j = 1:n_out
        switch outlier_side
            case -1
                [sample,rm_idx(j)] = min(temp);
                R(j) = (mean(temp,'omitnan')-sample);
            case 0
                [R(j),rm_idx(j)] = max(abs(temp-mean(temp,'omitnan')));
            case 1
                [sample,rm_idx(j)] = max(temp);
                R(j) = (sample-mean(temp,'omitnan')); 
        end
        R(j) = R(j)/std(temp,'omitnan');
        temp(rm_idx(j)) = NaN;
        
        p = 1-alpha/(n-j+1); 
        t = tinv(p,n-j-1);       
        lambda(j) = ((n-j)*t)/(sqrt((n-j-1+t^2)*(n-j+1)));
    end
    
    % And return a logical array of outliers
    idx = zeros(1,n);
    idx(rm_idx(1:(find(R>lambda,1,'last'))))=NaN;
    idx = ~isfinite(idx);
    
    if nargout > 1
        x2 = x(~idx);
    end
