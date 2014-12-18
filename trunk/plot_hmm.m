function [fig_handles fig_names fig_titles]=plot_hmm( S )

% [fig_handles fig_names fig_titles]=plot_hmm( S )
%
% Plots HMM diagnostics
%
% Inputs:
% hmm=S.hmm;
% block=S.block;
% tres=S.tres;
% Apca=S.Apca;
% NK=S.NK;
%
% MWW 2013

block=S.block;
tres=S.tres;
NK=S.NK;

fignum=1;

fig_handles(fignum)=figure; 
fig_names{fignum}='hmm_state_courses';
fig_titles{fignum}='HMM state time courses';
yticklab=[];
yticks=[];
for ii=1:NK,    
    ts=tres:tres:length(block(1).q_star)*tres;
    %ts=tres:tres:size(xin,2)*tres;        
    %ts=ts(delta*(M-1)+1:end);
    plot(ts,(block(1).q_star==ii)+2*(ii) ,'b');
    hold on;
    yticks=[yticks ii*2+0.5];
    yticklab{ii}=num2str(ii);
end;
set(gca,'YTick', yticks);
set(gca,'YTickLabel',yticklab);
plot4paper('time(s)','');
fignum=fignum+1;
    
if isfield(S,'hmm') && isfield(S,'Apca')
    Apca=S.Apca;
    hmm=S.hmm;

    fig_handles(fignum)=sfigure; 
    fig_names{fignum}='hmm_state_covs';
    fig_titles{fignum}='HMM state data covs';

    for ii=1:NK,

        cv=Apca*hmm.state(ii).Cov*Apca';
        snugplot(ceil(NK/2),2,ii); imagesc(abs(cv)); colorbar; 
    %        subplot(NK,1,(ii-1)*2+2); imagesc(abs(hmm.state(ii).Mu), [0 1]); colorbar; 

    end;
    fignum=fignum+1;
    
end;


StatePath = block(1).q_star;
%StatePath = simdata.Xclass(1:Nto);
NumberOccurences = zeros(1,NK);
FractionalOccupancy = zeros(1,NK);
MeanLifeTime = zeros(1,NK);

for ctr = 1:NK
    temp = StatePath == ctr;
    NumberOccurences(ctr) = sum(diff(temp) == 1);
    FractionalOccupancy(ctr) = sum(temp)/length(temp); 
    MeanLifeTime(ctr) = sum(temp)/NumberOccurences(ctr);
end
x=1:NK;

% plots
fig_handles(fignum)=sfigure; 
fig_names{fignum}='hmm_state_stats';
fig_titles{fignum}='HMM state stats';

subplot(1,2,1);bar(x,FractionalOccupancy*range(ts));
plot4paper('State #','Occupancy (s)');

subplot(1,2,2);bar(x,(tres)*MeanLifeTime);
plot4paper('State #','Mean life time (s)');



end

