function [ res ] = teh_graph_gmm_fit( S )

% [ res ] = teh_graph_gmm_fit( S )


try tmp=S.do_plots; catch, S.do_plots=1; end
try tmp=S.pvalue_th; catch, S.pvalue_th=0.05/length(S.data); end

res=[];

data=squash(S.data);

if S.do_fischer_xform
    logistic_data=log((data)./(1-data));
else
    logistic_data=data;
end

normalised_logistic_data=normalise(logistic_data);   

mix = gmm(1, 2, 'diag');
mix = gmminit(mix, normalised_logistic_data, foptions);
mix = gmmem(mix, normalised_logistic_data, foptions);

%%%%%%%%%%%%%%%%%%%%%
%% compute thresh
xx=min(normalised_logistic_data):range(normalised_logistic_data)/50:max(normalised_logistic_data);
xx2=linspace(0.003, 0.014, 50);

[prob] = gmmactiv(mix, xx');

prob2=prob;
for cc=1:2,
    prob2(:,cc)=prob2(:,cc)*mix.priors(cc); 
end

%if prob2(nearest(xx,0),1)<prob2(nearest(xx,0),2)
if mix.centres(2)<mix.centres(1)
     ind_null=2;
else
    ind_null=1;
end

indsearch=[nearest(xx,0):length(xx)];
th=xx(indsearch(nearest(prob2(indsearch,ind_null),S.pvalue_th)));

% check to see if 2nd dist is irrelevant if it is too similar to
% 1st dist
if abs(range(mix.centres))<1 || th<0
    th=[];
end

if isempty(th)
    th=max(normalised_logistic_data)+1;
end

orig_th=data(nearest(normalised_logistic_data,th));

%%%%%%%%%%%%%%%%%%%%%
%% do plots
if S.do_plots
    %%
    figure;
    set(gcf,'Position',[440   378   1000   420])

    subplot(131);
    plot(xx,prob2,'LineWidth',3)
    [prob3] = gmmprob(mix, xx');
    plot4paper('normalised edge strength','');

    subplot(132);
    hh=hist(normalised_logistic_data,xx);
    bar(xx,hh/sum(hh));
    ho;
    plot(xx,prob3/sum(prob3),'r','LineWidth',3);
    title(num2str(th));
    plot4paper('normalised edge strength','');
    
    subplot(133);
    [hh]=hist(data,xx2);
    bar(xx2,hh/sum(hh));
    plot4paper('coh','');
    a=axis;
    a(1)=xx2(1);   
    a(2)=xx2(end);
    axis(a)
    title(num2str(orig_th));

end

res.mix=mix;
res.normalised_th=th;
res.orig_th=orig_th;
res.data=normalised_logistic_data;

end

