function results = pacf(xin, min, max)

% results = pacf(xin, min, max)
%
% fits an AR(p) model to data xin. The model order p is
% automatically determined using the partial autocorrelation
% function method.
% min is minimum AR order, max is maximum AR order

N = size(xin);
results.order = -1;
x = xin;

for i = min:max
	y = x(i+1:N);% - mean(x);
	X = zeros(N(1)-i,i);
	
	for j = 1:i,
		X(:,j) = x((i+1-j):(N-j));% - mean(x);
	end;
    
    results.beta = y\X;
	if (abs(results.beta(i)) <  (1/N(1)) + 2/(N(1)^0.5) & results.order==-1 | i == max);
		results.order = i;
		break;
	end;
end;
