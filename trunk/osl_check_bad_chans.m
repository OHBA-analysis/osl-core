 function  [res fig_names fig_handles fig_titles]=osl_check_bad_chans(S)
   
% [res fig_names fig_handles]=osl_check_bad_chans(S)
%
% Used by osl_run_opt
%
% MWW 2013

res=[];
res.plot_fnames=[];

D=S.D;

if ~isfield(S, 'plot_basename'),
    S.plot_basename='';
end;
if ~isfield(S, 'plot_basetitle'),
    S.plot_basetitle=S.plot_basename;
end;

fig_titles={};
fig_names={};
fig_handles=[];
fig_index=1;

%%%%%%%%%%%%%%%%%%%%%%%%%

% all good chans in orig
chan_list=find(strcmp(chantype(D),S.modality));
chanind = chan_list;
Dx=D;

chanindDx = setdiff(chan_list, Dx.badchannels);

[c,ia,ib] = intersect(chan_list, Dx.badchannels);

badchans_in_modality=chan_list(ia);

if ~isempty(badchans_in_modality)
    disp(['Labelled bad ' S.modality ' chans=' mat2str(Dx.badchannels)]);
end;
    
chanindDx_inchanind=zeros(length(chanindDx),1);
for jj=1:length(chanindDx),
    tmp=find(chanind==chanindDx(jj));
    chanindDx_inchanind(jj)=tmp;            
end;

% all good trials in new
trials = Dx.indtrial(Dx.condlist,'good');

dat=D(chanind,:,trials);
datDx=Dx(chanindDx,:,trials);

datchan=feval('std',reshape(dat,size(dat,1),size(dat,2)*size(dat,3)),[],2);        
datchanDx=nan(size(datchan));
datchanDx(chanindDx_inchanind)=feval('std',reshape(datDx,size(datDx,1),size(datDx,2)*size(datDx,3)),[],2);

fig_handles(fig_index)=sfigure; 
fig_names{fig_index}=[S.plot_basename '_badchans_' S.modality];
fig_titles{fig_index}=[S.plot_basetitle 'badchans_' S.modality];
fig_index=fig_index+1;

snugplot(4,1,1);
[hs hsx]=hist(datchan,length(datchan)/10);
bar(hsx,hs);
title(['std' '(chans), ' S.modality]);
a1=axis;
snugplot(4,1,2);
plot(datchan,1:length(chanind),'*r');ho;
plot(datchanDx,1:length(chanind),'og');

a2=axis;
axis([a1(1) a1(2) a2(3) a2(4) ]);

snugplot(4,1,3);
[hs hsx]=hist(datchanDx,length(datchanDx)/10);
bar(hsx,hs);
a1=axis;
snugplot(4,1,4);
plot(datchanDx,1:length(chanind),'*g');
a2=axis;
axis([a1(1) a1(2) a2(3) a2(4) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_bch=badchans_in_modality;

if(length(list_bch)>0),
  
    % plot each time course for the badchans
    nbs=length(list_bch);
    
    if(nbs>16),
        
        jj=1;
        kk=1;
        fig_handles(fig_index)=sfigure; 
        fig_names{fig_index}=[S.plot_basename '_badchan' num2str(kk) '_' S.modality];
        fig_titles{fig_index}=[S.plot_basetitle 'Bad Chan' num2str(kk) ', ' S.modality];
        fig_index=fig_index+1;
                
        for ii=1:nbs,

            subplot(4,4,jj);
            data=D(list_bch(ii),:,1);
            plot(data(1:10:end));title(['chan ' num2str(list_bch(ii))]);
            
            if(jj==16 && ii~=nbs),

                jj=1;
                kk=kk+1;
                
                fig_handles(fig_index)=sfigure; 
                fig_names{fig_index}=[S.plot_basename '_badchan' num2str(kk) '_' S.modality];
                fig_titles{fig_index}=[S.plot_basetitle 'Bad Chan' num2str(kk) ', ' S.modality];
                fig_index=fig_index+1; 
                
            else
                jj=jj+1;
            end;
            
        end;
                        
    else,
        
        kk=1;
        fig_handles(fig_index)=sfigure; 
        fig_names{fig_index}=[S.plot_basename '_badchan' num2str(kk) '_' S.modality];
        fig_titles{fig_index}=[S.plot_basetitle '_badchan' num2str(kk) '_' S.modality];
        fig_index=fig_index+1;
        for jj=1:nbs,
            subplot(ceil(sqrt(nbs)),ceil(sqrt(nbs)),jj);
            data=D(list_bch(jj),:,1);
            plot(data(1:10:end));title(['chan ' num2str(list_bch(jj))]);
        end;


    end;
end;
