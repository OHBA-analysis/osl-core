function h = spectrograms(D)


    for mm=1:length(S.modalities),
        fig_handles=sfigure;
        set(fig_handles,'Position',[1 1 1300 450]);

        S=[]; 
        S.D=D_continuous; 
        S.chantype=opt.modalities{mm};
        S.do_plot=false;
        S.cut_badsegments=true;
        [spect,F,T]=osl_plotspectrogram(S);

        %do plot
        imagesc(T,F,spect); colorbar
        set(gca,'ydir','normal');
        plot4paper('time (s)','freq (hz)');

        report=osl_report_set_figs(report,['spectogram_' opt.modalities{mm}],fig_handles,['DATA: spectogram for ' opt.modalities{mm}]);
        report=osl_report_print_figs(report);
    end



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

    if isobject(Sin)
        D=Sin;
        Sin=[];
        Sin.D=D;
    end

    try D=Sin.D; catch error('Must specify S.D'); end
    try Sin.do_plot=Sin.do_plot; catch Sin.do_plot=true; end
    try D.fsample; catch D=spm_eeg_load(D); end
    try chantype=Sin.chantype; catch chantype='MEGGRAD'; end
    try chaninds=Sin.chaninds; catch chaninds=1:D.nchannels; end
    try cut_badsegments=Sin.cut_badsegments; catch cut_badsegments=0; end

        
    if D.ntrials>1
        error('Only works on continuous data')
    end

    goodsamples = find(good_samples(D));

    if cut_badsegments
        [S,F,T] = plotspectrogram(D(1,goodsamples,1),512,512*0.75,1024,D.fsample);
    else
        [S,F,T] = plotspectrogram(D(1,:,:),512,512*0.75,1024,D.fsample);
    end

    chindex = 1:D.nchannels;
    ch = chindex(find(strcmp(D.chantype,chantype)));
    %ch=ch(chaninds);
    %ch=chaninds(ch);

    S(:)=0;
    for i = 1:length(ch)
      if cut_badsegments
        S=S+plotspectrogram(D(ch(i),goodsamples,1),512,512*0.75,1024,D.fsample);
      else
        S=S+plotspectrogram(D(ch(i),:,:),512,512*0.75,1024,D.fsample);
      end
      
      %P{i} = mean(S{i},2);
      %disp([num2str(i) '/' num2str(length(ch))])
    end
    S=S/length(ch);

    %% Plot mean spectrogram
    if Sin.do_plot    
        imagesc(T,F,S); colorbar
        set(gca,'ydir','normal');
        plot4paper('time (s)','frequency (Hz)');
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


function [varargout] = osl_plotspectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)

    [S,F,T,P] = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs);
    S = log10(P);


    switch nargout,
     case 0,
       imagesc(T,F,S);
       set(gca,'ydir','normal');
       xlabel('time (s)');
       ylabel('frequency (Hz)');
     case 1
       varargout = {S};
     case 2
       varargout = {S,F};
     case 3
       varargout = {S,F,T};
    end

    end
