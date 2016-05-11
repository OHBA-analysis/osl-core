function h = red2yellow(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end

r=ones(m,1);
g=[0:1/(m-1):1]';
b=zeros(m,1);

h=[r g b];
