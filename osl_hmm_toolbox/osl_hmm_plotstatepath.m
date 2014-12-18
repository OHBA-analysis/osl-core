function osl_hmm_plotstatepath(hmm)
% Plots the HMM statepath
% osl_hmm_plotstatepath(hmm)
% AB 2013

if ~isfield(hmm,'fsample') || isempty(hmm.fsample)
    fs = 1;
else
    fs = hmm.fsample;
end

cla(gca); hold(gca,'on')

k = max(hmm.statepath);

col = colormap(gca,'lines');

for s = 1:k
  plot((1:length(hmm.statepath))/fs, 0.8*double(hmm.statepath==s)+ s - 0.4,'color',col(s,:));
end

if ~isfield(hmm,'fsample') || isempty(hmm.fsample)
  xlabel('Samples');    
else
  xlabel('Time (s)');
end
 
ylabel('State #')
 
set(gca,'xlim',[0 length(hmm.statepath)/fs])
set(gca,'ylim',[0 hmm.K]+0.5)


end
