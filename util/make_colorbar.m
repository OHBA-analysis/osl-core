function make_colorbar( cs, cmap, dec_points )

% make_colorbar( cs, cmap )
%
% e.g. cs=[20,26];cmap=autumn(64);
%
% MWW 2013

if nargin<2,
    cmap=autumn(64);
end;

if nargin<3,
    dec_points=1;
end;

num=4;
cs=cs(1):(cs(end)-cs(1))/(num-1):cs(end);
cs=round(cs*10^dec_points)/10^dec_points; % just show one decimal point

imagesc(flipud(repmat([1:64],7,1)'));
set(gca,'Xtick',[]);
set(gca,'Ytick',[1:63/(num-1):64]');
set(gca,'Yticklabel',fliplr(cs));
set(gca,'FontSize',14);

colormap(gca,cmap);
axis equal; axis tight;
