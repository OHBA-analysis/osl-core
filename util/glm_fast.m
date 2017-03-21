function [tstat, pe, dof, cope, varcope]=glm_fast(y,x,pinvxtx,pinvx,cs,con,abs_pes)

% y = demean(randn(100,1))+s'*0.1+s2'*2;

% t=[1:100];s=demean(sin(2*pi*t/20));s2=demean(sin(2*pi*t/30));x=[s; s2]'; c=[1 1]';x=demean(x);[tstat, pe]=glm(y,x,c)

% t=[1:100];s=demean(sin(2*pi*t/20));s2=demean(sin(2*pi*t/30));x=[s+s2]'; c=[1]';x=demean(x);[tstat, pe]=glm(y,x,c)

if(size(y,1)~=size(x,1))
    error('size(y,1)~=size(x,1)');
end;

pe=pinvx*y;

if abs_pes,
    pe=abs(pe);
end;

r=y-x*pe;

vr=diag(r'*r/(size(y,1)-size(x,2)));

if(con>0),
  c=cs{con}';

  varcope=c'*pinvxtx*c*vr;
  
  cope=c'*pe;

  tstat=cope./sqrt(varcope)';
else
    for cc=1:length(cs),
        c=cs{cc};
 
        varcope(cc)=c'*pinvxtx*c*vr;
        cope(cc)=c'*pe;
        tstat(cc)=cope(cc)/sqrt(varcope(cc));  
    end;
end;

dof=length(y)-size(x,2);