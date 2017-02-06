function [cope,varcope,tstat]=ols(data,des,tc)
% [COPE,VARCOPE]=OLS(DATA,DES,TC)
% DATA IS T x V
% DES IS T x EV (design matrix)
% TC IS NCONTRASTS x EV  (contrast matrix)

if(size(data,1)~=size(des,1))
  error('OLS::DATA and DES have different number of time points');
elseif(size(des,2)~=size(tc,2))
  error('OLS:: DES and EV have different number of evs')
end


pdes=pinv(des);
prevar=diag(tc*pdes*pdes'*tc');
R=eye(size(des,1))-des*pdes;
tR=trace(R);
pe=pdes*data;
cope=tc*pe; 
if(nargout>1)
  res=data-des*pe;
  sigsq=sum(res.*res/tR);       
  varcope=prevar*sigsq;
  if(nargout>2)
    tstat=cope./sqrt(varcope);
  end
end