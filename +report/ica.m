function h = ica(D)
    % Returns handles to figures summarizing ICA results

    % Couple of different ICA type plots

    % We will plot - any ICs that were bad, and any which were at one point automatically marked bad
    % Possible states are combinations of - auto=true/false, manual=true/false
    if ~isfield(D,'ica')
        error('MEEG object is missing ICA results - make sure osl_africa has been run first');
    end
    
    h = [];

    auto_marked = find(~cellfun(@isempty,D.ica.auto_reason));
    to_plot = union(D.ica.bad_components,auto_marked);

    % Reconstruct timecourses
    D = D.montage('switch',0);
    samples_of_interest = good_samples(D);
    samples_of_interest = reshape(samples_of_interest,1,D.nsamples*D.ntrials);
    sm = D.ica.sm(D.ica.chan_inds,:);
    tc = (D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))').';
    tc(:,~samples_of_interest) = NaN;
    t = (1:size(tc,2))./D.fsample;

    % If topos not precomputed, then compute them now
    if isempty(D.ica.topos)
        fprintf('No precomputed topos present, computing them now...\n')
        topos = [];
        modalities = D.ica.modalities; 
        for m = 1:numel(modalities)
            topos = [topos component_topoplot(D,D.ica.sm,modalities(m))];
        end
    else
        topos = D.ica.topos;
    end

    % PLOT EIGENSPECTRA
    h(end+1) = figure('name','ICA Eigenspectra','tag','ica_eigenspectra');
    semilogy(D.ica.eigs_preNorm);
    hold on
    semilogy(D.ica.eigs_postNorm,'r--');
    title('Raw and normalised eigen spectra'); legend('Raw', 'Normalised');

    % PLOT BAD COMPONENTS
    n_cols = length(D.ica.modalities);
    for j = 1:length(to_plot)
        h(end+1) = figure('name',sprintf('IC component %d',to_plot(j)),'tag',sprintf('ic_component_%d',to_plot(j)));
        pos = get(gcf,'Position');
        set(gcf,'Position',pos.*[1 1 1 1.5])

        ax = [];
        for k = 1:n_cols
            ax(k) = subplot(3,n_cols,k);;
            struct2handle(topos(to_plot(j),k).children,ax(k)); ;
            set(gca,'CLim',[-1 1]*max(abs(get(gca,'CLim'))));;
            colormap(osl_colormap('rwb'));
            axis equal
            axis tight
            axis off
            colorbar
        end

        ax(n_cols+1) = subplot(3,2,[n_cols+1:n_cols*2]);
        l = plot(t,tc(to_plot(j),:));
        if any(to_plot(j)==D.ica.bad_components) % If the component was actually removed
            set(l,'Color',[204     0     0] / 255);
            title(sprintf('IC%d - REMOVED',to_plot(j)))
        else
            set(l,'Color',[ 48   128    20] / 255);
            title(sprintf('IC%d - NOT REMOVED',to_plot(j)))
        end
        xlabel('Time (s)')
        ylabel('Signal')
        set(gca,'XLim',[min(t) max(t)]);

        ax(n_cols+1) = subplot(3,2,[n_cols*2+1:n_cols*3]);
        [P,f] = pwelch(tc(to_plot(j),~isnan(tc(to_plot(j),:))),1024,512,1024,D.fsample);
        l = semilogy(f,P,'Color',[ 48   128    20] / 255);
        if any(to_plot(j)==D.ica.bad_components) % If the component was actually removed
            set(l,'Color',[204     0     0] / 255);
        end
        set(gca,'XLim',[min(f) max(f)]);
        xlabel('Frequency (Hz)')
        ylabel('Power spectral density');

        if isempty(D.ica.auto_reason)
            title('Auto rejection not performed');
        elseif any(to_plot(j)==auto_marked)
            title(sprintf('Automatically marked bad:%s',D.ica.auto_reason{to_plot(to_plot(j)==auto_marked)}),'interpreter','none');
        else
            title(sprintf('Not automatically marked bad'));
        end

    end

    % PLOT REJECTION VS METRIC
    fields = fieldnames(D.ica.metrics);
    is_rejected = ismember(1:size(tc,1),D.ica.bad_components);

    for j = 1:length(fields)

        data = D.ica.metrics.(fields{j}).value;
        dataclean=data(~is_rejected);
        [~,ia] = sort(data,'descend');
        ib(ia)=1:length(data);
        dataclean = sort(dataclean,'descend');

        h(end+1) = figure('name',sprintf('ICA metric: %s',fields{j}),'tag',sprintf('ica_metric_%s',fields{j}));

        subplot(2,2,1)
        [hs hsx]=hist(data,length(data)/10);        
        bar(hsx,hs);
        xlabel(fields{j},'interpreter','none');
        box on

        subplot(2,2,2)
        colormap(gca,[ 48   128    20;204     0     0] / 255);
        scatter(ib,data,50,'*','CData',1+is_rejected);
        xlabel(fields{j},'interpreter','none');
        ylabel('Ordered IC')
        box on

        subplot(2,2,3)
        [hs hsx]=hist(dataclean,length(dataclean)/10);        
        bar(hsx,hs);
        xlabel(fields{j},'interpreter','none');
        box on

        subplot(2,2,4)
        scatter(dataclean,1:length(dataclean),'*','CData',[ 48   128    20] / 255);
        xlabel(fields{j},'interpreter','none');
        ylabel('Ordered IC')
        box on
    end



