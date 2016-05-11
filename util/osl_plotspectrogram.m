function [S,F,T] = osl_plotspectrogram(Sin)

% [spectrogram,freq,time] = osl_plotspectrogram(S)
%
% Inputs:
% S.D
% S.chantype (e.g. 'MEGGRAD','MEGPLANAR')
% S.do_plot
%
% e.g. 
% S=struct('D',D);
% [spectrogram,F,T] = osl_plotspectrogram(S); 
% figure;imagesc(T,F,spectrogram); colorbar
%
% AB & MWW 2016

try D=Sin.D; catch error('Must specify S.D'); end
try Sin.do_plot=Sin.do_plot; catch Sin.do_plot=true; end
try D.fsample; catch D=spm_eeg_load(D); end

try chantype=Sin.chantype; catch chantype='MEGGRAD'; end

[S,F,T] = plotspectrogram(D(1,:,:),512,512*0.75,1024,D.fsample);

chindex = 1:D.nchannels;
ch = chindex(strcmp(D.chantype,chantype));

S(:)=0;
for i = 1:length(ch)
  S=S+plotspectrogram(D(ch(i),:,:),512,512*0.75,1024,D.fsample);
  %P{i} = mean(S{i},2);
  %disp([num2str(i) '/' num2str(length(ch))])
end
S=S/length(ch);

%% Plot mean spectrogram
if Sin.do_plot
    figure
    imagesc(T,F,S,[0 3]); colorbar
    set(gca,'ydir','normal');
    xlabel('time (s)');
    ylabel('frequency (Hz)');
end

if 0
    
    % Get scalp layout
    chan_ind = D.selectchannels(chantype);
    p = D.coor2D(chan_ind);

    %% Plot data in x & y on scalp layout
    f1 = figure;

    sens_offset = zeros(size(chan_ind));
    sens_offset(1:2:end) = -0.025;
    sens_offset(2:2:end) = 0;

    clear ax
    for n = 1:length(ch)
      ax(n) = axes('position',[p(1,n)-sens_offset(n),p(2,n)-0.025,0.025,0.05]);
      plot(F,mean(P{n},2))
      set(gca,'ydir','normal');
      set(ax(n),'xtick',[]);
      set(ax(n),'ytick',[]);
    end


    %% Run this to plot select individual sensors to plot separately
    figure(f1)
    [~,~,but] = ginput(1);
    ind = find(ax==gca);
    figure
    imagesc(T,F,S{ind});
    title(D.chanlabels(ind))
    set(gca,'ydir','normal');
    xlabel('time (s)');
    ylabel('frequency (Hz)');
end
