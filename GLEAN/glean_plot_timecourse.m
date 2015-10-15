function glean_plot_timecourse(GLEAN)
% Plots GLEAN inferred time courses (state time courses for HMM and
% component time courses for ICA)

load(GLEAN.model.model)

model = char(fieldnames(GLEAN.model.settings));
switch model
    case 'hmm'
        load(GLEAN.model.model);
        tc = cell2mat(arrayfun(@(k) hmm.statepath == k,1:hmm.K,'UniformOutput',0))';
        plotSpacing = 1.2*ones(size(tc,1),1);
        plotCentres = 0.5*ones(size(tc,1),1);
        ylbl = 'State';
    case 'ica'
        load(GLEAN.model.model);
        tc = ica.tICs;
        plotSpacing = 6*std(tc,[],2);
        plotCentres = mean(tc,2);
        ylbl = 'Component';
end

% Plot multiple timeseries
figure('color','w'), hold on, box on
plot((tc + repmat(cumsum(plotSpacing),1,size(tc,2)))')
axis tight
set(gca,'ylim',get(gca,'ylim') + ([-1 1] .* 0.01*diff(get(gca,'ylim'))));
sessionMarkers = find(diff(subIndx));
stem(sessionMarkers,max(get(gca,'ylim'))*ones(numel(GLEAN.data)-1,1),'--k','marker','none')

xlabel('Session','fontsize',14)
xtick = [0 sessionMarkers] + 0.5*diff([0 sessionMarkers size(tc,2)]);
set(gca,'xtick',xtick,'xticklabel',1:length(sessionMarkers)+1)

ylabel(ylbl,'fontsize',14)
ytick = cumsum(plotSpacing) + plotCentres;
set(gca,'ytick',ytick,'yticklabel',1:size(tc,1))
set(gca,'ticklength',[0 0])

end