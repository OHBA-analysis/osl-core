function [fig_handles] = oat_plot_hmm_states( oat )

% fig_handles=oat_plot_hmm_states( oat )

if strcmp(oat.source_recon.modalities{1},'EEG')
    modality_meeg='EEG';
else
    modality_meeg='MEG';
end

fig_handles=[];

%%%%%%%%%%%%%%%%%%%%%%%
%% plot hmm classes

res=oat_load_results(oat,oat.source_recon.results_fnames{1});

D=oat_get_sensordata(res);

block=res.block;
samples2use=res.samples2use;

classchanind=find(strcmp(D.chanlabels,'Class'));
hmm_class=D(classchanind,find(res.samples2use),res.trials);
hmm_class2=reshape(hmm_class,size(hmm_class,1),size(hmm_class,2)*size(hmm_class,3));

NK=max(unique(hmm_class))

if D.ntrials==1
    tim=D.time(find(samples2use));
else
    tres=1/D.fsample;
    tim=tres:tres:tres*length(hmm_class2);
end;

if(isfield(res,'pca_tra'))
    chans=find(strcmp(D.chantype,'MEGPCACOMP'));
    data=((res.pca_tra'*reshape(D(chans,find(res.samples2use),res.trials),length(chans),(numel(find(res.samples2use))*length(res.trials))))');
else
    chans=find(strcmp(D.chantype,modality_meeg));
    data=((reshape(D(chans,find(res.samples2use),res.trials),length(chans),(numel(find(res.samples2use))*length(res.trials))))');
end;

data2=reshape(data',size(data,2),numel(find(res.samples2use)),length(res.trials));
data2=mean(abs(data2),1);

for cond=1:length(oat.source_recon.conditions),
    trls{cond}=intersect(D.indtrial(oat.source_recon.conditions{cond},'good'),res.trials);
    hmm_class_cond{cond}=D(classchanind,find(res.samples2use),trls{cond});
end;
    
    
for ii=1:NK,

    tmp=zeros(length(oat.source_recon.conditions),size(hmm_class,2));

    for cond=1:length(oat.source_recon.conditions),
        tmp(cond,:)=sum(hmm_class_cond{cond}==ii,3)/size(hmm_class_cond{cond},3);
    end;
    fig_handles{length(fig_handles)+1}=figure;
    plot(D.time(res.samples2use),tmp,'LineWidth',2);plot4paper('time (s)','frac occupancy');   
    title(['State #' num2str(ii)]); 
    legend(oat.source_recon.conditions);
    
end;

tmp=zeros(NK,size(hmm_class,2));
tmpdat2=zeros(NK,size(hmm_class,2),size(hmm_class,3));

fig_handles{1}=figure; hold on;
for ii=1:NK,
    
    tmp(ii,:)=sum(hmm_class==ii,3)/size(hmm_class,3);
    
    tmpdat=data2;
    tmpdat(hmm_class~=ii)=0;
    tmpdat2(ii,:,:)=tmpdat;
    
    plot(tim,(hmm_class2==ii)+2*(ii) ,'b');
    ylabs{ii}=num2str(ii);

end;
set(gca,'YTick',[2.5:2:2.5+NK*2]);
set(gca,'YTickLabel',ylabs)
plot4paper('time (secs)','State #');

fig_handles{length(fig_handles)+1}=figure;plot(D.time(res.samples2use),tmp,'LineWidth',2);
plot4paper('time (s)','frac occupancy');

%fig_handles{length(fig_handles)+1}=figure;plot(D.time(res.samples2use),mean(tmpdat2,3),'LineWidth',2);plot4paper('time (s)','HMM class no.');
%title('Av data');

fig_handles{length(fig_handles)+1}=figure;
for ii=1:NK,
    subplot(NK,1,ii);
    plot(D.time(res.samples2use),squeeze(tmpdat2(ii,:,:)),'LineWidth',1);plot4paper('time (s)','HMM class no.');
end;

fig_handles{length(fig_handles)+1}=figure;ho;

if(isfield(res,'pca_tra'))
    chans=find(strcmp(D.chantype,'MEGPCACOMP'));
    data=((res.pca_tra'*reshape(D(chans,find(res.samples2use),res.trials),length(chans),(numel(find(res.samples2use))*length(res.trials))))');
else
    chans=find(strcmp(D.chantype,modality_meeg));
    data=((reshape(D(chans,find(res.samples2use),res.trials),length(chans),(numel(find(res.samples2use))*length(res.trials))))');
end;

data=mean(abs(data),2);
%data=data(:,:);
sqdata=squash(data);
rng=max(sqdata)-min(sqdata);

ii=0;
tmp=data;
plot(tim,tmp+rng*(ii) ,'b');
for ii=1:NK,

    tmp=data;
    tmp(hmm_class2~=ii,:)=0;
    plot(tim,tmp+rng*(ii) ,'b');

end;

plot4paper('time (secs)','States');
    
end

