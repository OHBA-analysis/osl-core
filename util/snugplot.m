function snugplot(a,b,c,pc)
% snugplot(a,b,c,pc)
% like subplot(a,b,c)
% pc is the fraction of white space along each axis.
% TB04

if(c>(a*b))
  error('ya twat');
end
if(nargin==3)
  pc=0.25;
end


numplots=a*b;
percy=(1-pc)/a;
percx=(1-pc)/b;
ygap=pc/(a+1);
xgap=pc/(b+1);


row=ceil(c/b);
row=a-row+1;
col=c-(ceil(c/b)-1)*b;

subplot('position',[xgap+(xgap+percx)*(col-1) ygap+(ygap+percy)*(row-1) percx percy]);

set(gcf,'PaperPositionMode','auto');

