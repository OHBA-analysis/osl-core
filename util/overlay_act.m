function actminmax=overlay_act(act, bg, cmapname, cbar, actminmax, bgminmax, deact_cmapname, deactminmax)

% overlay_act(act, bg, cmapname, cbar, actminmax, bgminmax, deact_cmapname, deactminmax)

deact=-act;

if nargin<3
  cmapname='red2yellow';
end;

if nargin<8
  deactminmax=[max(squash(deact))*2 max(squash(deact))*3]; % do not show any -ve
end;

if nargin<7
  deact_cmapname='bone';
end;

if nargin<6 || isempty(bgminmax)
  bgminmax=[min(squash(bg)),max(squash(bg))];
end;

if nargin<5
  actminmax=[max(0,min(squash(act))), max(squash(act))];
end;

if nargin<4
  cbar=0;
end;

%if(actminmax(1)<min(squash(act))), actminmax(1)=min(squash(act)); end;
%if(actminmax(2)>max(squash(act))), actminmax(2)=max(squash(act)); end;

%if(deactminmax(1)<min(squash(deact))), deactminmax(1)=min(squash(deact)); end;
%if(deactminmax(2)>max(squash(deact))), deactminmax(2)=max(squash(deact)); end;

exec_str=strcat('[gray(100);',cmapname,'(100);',deact_cmapname,'(100)];');
cmap=eval(exec_str);

bg2=(bg-bgminmax(1))./range(bgminmax);
bg2(bg2<0)=0;
bg2(bg2>1)=0.999;

%bg2=[0.5 0.5 0.5 0.5 0.5];act=[-20 -10 10 50 300];deact=-act;

act2=(act-actminmax(1))./range(actminmax)+1;
act2(act2<1)=0.999;
act2(act2>=2)=1.999;
act2(act<0)=0.999;

deact2=(deact-deactminmax(1))./range(deactminmax)+2;
deact2(deact2<2)=0.999;
deact2(deact2>=3)=3;
deact2(deact<0)=0.999;

bg2(act2>1)=act2(act2>1);

bg2(deact2>2)=deact2(deact2>2);

%snugplot(1,1,1,0.02);

imagesc(bg2,[0 3]);
axis image;
axis off;
colormap(gca,cmap);

%tmp=get(gcf,'Position');
%tmp(3)=tmp(4);
%set(gcf,'Position',tmp);

if(cbar>0),
  fh=gcf;
  figure;
  tmp=actminmax(1):range(actminmax)/100:actminmax(2);
  h=imagesc(repmat(tmp,10,1)');
  axis xy;
  axis image;
  colormap(cmapname);
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  text(2,-3,sprintf('%1.1f',actminmax(1)),'FontSize',12);
  text(2,105,sprintf('%1.1f',actminmax(2)),'FontSize',12);
  tmp=get(gcf,'Position');
  tmp(3)=50;
  set(gcf,'Position',tmp);
  set(gcf,'PaperPositionMode','auto');
  figure(fh);
end;
