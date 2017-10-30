% Plot if required

D.ica.modalities

if S.do_plots && ~isempty(bad_components)
    for j = 1:length(bad_components)
        figs.handles(end+1) = figure('Position',[1 1 1500 1000]);
        figs.names{end+1} =  sprintf('IC %d - %s',bad_components(j),artefacts{bad_components(j)});
        figs.titles{end+1} = sprintf('AFRICA: IC %d - %s',bad_components(j),artefacts{bad_components(j)});    

        c=1;
        ncols=max(numel(modalities),2);
                   
        for m =1:numel(modalities)
            subplot(2,ncols,c); 
            component_topoplot(D,D.ica.sm(:,bad_components(j)),modalities{m},true);
            title(modalities{m}); 
            axis tight;
            %title(['IC' num2str(bad_components(j)) ', ' confirmed_artefact_names{ii} ' Artefacts' modalities(m)]);%axis tight;
            c=c+1;
        end
        
        c=1;
        
        subplot(2,ncols,ncols+c);
        plot(freq_ax,abs_ft(bad_components(j),1:floor(length(abs_ft)/2)));title('Frequency Spectrum');
        xlabel('Frequency (Hz)');
        ylabel('Spectrum');
        axis tight;
        c=c+1;

        subplot(2,ncols,ncols+c);

        if D.ntrials == 1
            t = D.time;
        else
            t = (1:length(samples_of_interest))./D.fsample;
        end
        
        plot(t,tc(bad_components(j),:));
        title(sprintf('Component %d - %s',bad_components(j),artefacts{bad_components(j)}));
        xlabel('Time (s)');
        axis tight;
            
    end
end




%% RANK AND PLOT

if do_plots
    fig_handles(cc)=fig_handle;
    fig_names{cc}=fig_name;
    fig_titles{cc}=fig_title;
end

cc=cc+1;

% if strcmp(direction,'descend'),
%     num_comps2reject=max(comps2reject);
% else
%     num_comps2reject=min(comps2reject);
% end
% 
% num_comps=length(comps2reject);
%         
% if(num_comps>max_num_artefact_comps),
%     warning(['AFRICA detected ' num2str(num_comps) ' ICs with ' figlab ' > ' num2str(thresh) '. Only rejecting the ' num2str(max_num_artefact_comps) ' worst components.']);
% 
%     num_comps=max_num_artefact_comps;
% end
% 
% comps2reject=comps2reject(1:num_comps);

if ~do_plots
    return
end

for i = 1:length(comps2reject),

    fig_handles(cc)=sfigure; set(fig_handles(cc),'Position',[1 1 1500 1000]);
    fig_names{length(fig_names)+1}=['artefact_high_kurtosis_ic_' num2str(comps2reject(i))];
    fig_titles{length(fig_titles)+1}=['AFRICA: high kurtosis IC no.' num2str(comps2reject(i))];
    
    cc=cc+1;

    ncols=max(numel(modalities),2);

    for m = 1:numel(modalities)
        subplot(2,ncols,m);
        component_topoplot(D,sm(:,comps2reject(i)),modalities{m},true);
        title({['Component ' num2str(comps2reject(i)) ': ' figlab ' = ' num2str(met_ord(i))] modalities{m}}); axis tight;
    end
    
    if D.ntrials == 1
        t = D.time;
    else
        t = (1:length(samples_of_interest))./D.fsample;
    end

    subplot(2,ncols,ncols+1); plot(t(samples_of_interest),tc(comps2reject(i),:)); title({['Component ' num2str(comps2reject(i)) ': ' figlab ' = ' num2str(met_ord(i))] 'Independent Time Course'}); xlabel('Time (s)');axis tight;
    subplot(2,ncols,ncols+2); plot(freq_ax,abs_ft(comps2reject(i),1:floor(length(abs_ft)/2))); title({['Component ' num2str(comps2reject(i)) ': ' figlab ' = ' num2str(met_ord(i))] 'Frequency Spectrum'}); xlabel('Frequency (Hz)');axis tight;

end



%% Find outliers
if(num_comps>max_num_outliers),
 warning(['AFRICA detected ' num2str(num_comps) ' ' fig_label ' ICs. Only rejecting the ' num2str(max_num_outliers) ' worst components.']);

 num_comps=max_num_outliers;
end

outlier_inds=outlier_inds_sorted(1:num_comps);

 if(do_plots)    
    fig_handle=sfigure;
    fig_name=['hist_' fig_label];
    fig_title=['AFRICA: histogram of ' fig_label];
    
    subplot(221);        
    [hs hsx]=hist(data,length(data)/10);        
    bar(hsx,hs);
    plot4paper(fig_label,'');
    a1=axis;
    subplot(222);
    plot(data,1:length(data),'*g');ho;
    plot(data(outlier_inds),outlier_inds,'or');
    plot4paper(fig_label,'ordered IC');
    hold off;
    a2=axis;
    %axis([a1(1) a1(2) a2(3) a2(4) ]);

    subplot(223);     
    dataclean=data;
    dataclean(outlier_inds)=[];
    [hs hsx]=hist(dataclean,length(dataclean)/10);        
    bar(hsx,hs);
    plot4paper(fig_label,'');
    a1=axis;
    subplot(224);
    plot(dataclean,1:length(dataclean),'*g');ho;
    plot4paper(fig_label,'ordered IC');
    hold off;
    a2=axis;
end