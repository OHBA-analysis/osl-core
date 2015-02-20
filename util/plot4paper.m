function plot4paper(xtext,ytext,fsize)

if nargin<3
    fsize=16;
end;

a=axis;
if(nargin<2)ytext='';end;
set(gca,'fontsize',fsize);
set(gca,'LineWidth',2);
xlabel(xtext,'fontsize',fsize);
ylabel(ytext,'fontsize',fsize);
axis(a);