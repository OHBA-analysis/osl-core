% BOXCOX: Box-Cox maximum-likelihood transformation to find the optimal power 
%         transformation of an empirical distribution to symmetry.  The transformed 
%         data are scaled to have the same mean and standard deviation as the original 
%         data.  Ignores missing data.
%
%     Syntax:  [xp,lambda,c,Lmax] = boxcox(x,{plotflag})
%
%         x =        vector of observations for a single variable.
%         plotflag = optional boolean flag producing, if true, histograms of
%                      original and transformed data and plot of L vs lambda 
%                      [default = 0].
%         -------------------------------------------------------------------------
%         xp =       corresponding vector of transformed variable.
%         lambda =   Box-Cox parameter.
%         c =        value added to data before transforming to ensure all positive 
%                      values.
%         Lmax =     max log-likelihood function value.
%

% RE Strauss, 2/10/97
%    9/3/99 - changed plot colors for Matlab v5.
%   12/8/02 - change fmin() to fminbnd(); output the constant c;
%             added error message for missing data.
%   1/2/03 -  added optional histograms.
%   1/3/03 -  ignore missing data; produce separate histograms rather than subplots. 

function [xp,lambda,c,Lmax] = boxcox1(x,plotflag)
  if (nargin < 1) help boxcox; return; end;
  
  if (nargin < 2) plotflag = []; end;

  if (isempty(plotflag)) plotflag = 0; end;
  
  norig = length(x);
  i = find(isfinite(x));                  % Remove missing data
  x = x(i);
  n = length(x);

  xmin = min(x);
  if (xmin <= 0)                          % Distribution must be positive
    c = abs(xmin)+1;
    x = x+c;
  else
    c = 0;
  end;

  xmean = mean(x);
  xstd = std(x);

  lambda = fminbnd(@boxcoxf,-10,10,[],x);  % Optimize lambda
  
  if (abs(lambda) > eps)                
    xp = ((x.^lambda)-1)/lambda;
  else
    xp = log(x);
  end;

  Lmax = -((n-1)/2)*log(var(xp)) + (lambda-1)*((n-1)/n)*sum(log(x));

  if (plotflag)
    figure;
    histgram(x-c);
    putxlab('Data');
    puttitle('Original data');
    
    figure;
    histgram(xp-c);
    putxlab('Data');
    puttitle('Transformed data');
    
    figure;
    lvect = linspace(lambda-2,lambda+2);
    Lvect = lvect;
    for j=1:length(lvect)
      Lvect(j) = -boxcoxf(lvect(j),x);
    end;

    xmin = min(lvect) - 0.05*range(lvect);
    xmax = max(lvect) + 0.05*range(lvect);
    ymin = min(Lvect);
    ymax = max(Lvect) + 0.05*range(Lvect);

    plot(lvect,Lvect,'k');
    axis([xmin xmax ymin ymax]);
    putxlab('lambda');
    putylab('Log-likelihood (L)');
    hold on;
    plot([lambda lambda],[ymin Lmax],'k:');
    hold off;
  end;
  
  xpp = xp;
  xp = NaN*ones(norig,1);
  xp(i) = xpp;

  return;

end

% BOXCOXF: Objective function for boxcox().  x must be strictly positive.
%         This application minimizes L rather than maximizes -L.
%

% RE Strauss, 2/10/97
%   12/8/02 - replace nu by n-1.

function L = boxcoxf(lambda,x)
  n = length(x);

  if (abs(lambda) > eps)
    xp = ((x.^lambda)-1)/lambda;
  else
    xp = log(x);
  end;

  L = -((n-1)/2)*log(var(xp)) + (lambda-1)*((n-1)/n)*sum(log(x));
  L = -L;

  return;
end
