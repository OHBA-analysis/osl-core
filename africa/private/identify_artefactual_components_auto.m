%% osl_africa.m
%
% AfRICA - ArteFact Rejection using Independent Component Analysis
%
% Syntax: [fname_out scanner_artefacts_chan]=osl_africa(S)
% S needs to contain be produced inside OSL_AFRICA.m
% OPTIONAL inputs:
%   -  do_plots: set to 1 to output summary plots of artefact components.
%      set to 0 to switch off plotting.DEFAULT = 1.
%   -  ident.do_kurt: set to positive integer "n" to view and select artefacts based on 
%      "n" low/high kurtosis. Set to -1 to view all tCIs in desending
%      order.
%      DEFAULT = 0.
%   -  ident.do_mains: set to 1 to remove 50Hz mains artefacts.
%      set to 0 to leave mains components.DEFAULT = 1.
% Output:
%   - bad_components: A list of the components identified as bad.
%
% HL+MWW 2013

function [bad_components,metrics,figs] = identify_artefactual_components_auto(D,S)

    if isscalar(S.artefact_chans_corr_thresh)
        S.artefact_chans_corr_thresh = S.artefact_chans_corr_thresh*ones(size(S.artefact_channels));
    end
    
    % Compute metrics
    D = D.montage('switch',0);
    
    [metrics,tc] = compute_metrics(D,S.do_mains,S.mains_frequency,S.do_kurt,S.do_cardiac,S.artefact_channels);
    artefacts = cell(size(tc,1),1); % Record artefacts. Empty content means no artefact, otherwise string with reason for rejection

    modalities = unique(D.chantype(find(strncmpi(S.modality,D.chantype,3))));  % changed by DM

    %% Select only good data for classification
    samples_of_interest = ~all(badsamples(D,':',':',':'));
    samples_of_interest = reshape(samples_of_interest,1,D.nsamples*D.ntrials);

    num_ics = D.ica.params.num_ics;
    abs_ft  = abs(fft(demean(tc(:,samples_of_interest),2),[],2)); % Not entirely happy about using only samples_of_interest but this is what was used previously
    freq_ax = 0:(D.fsample/2)/(floor(length(abs_ft)/2)-1):(D.fsample/2);

    minhz = 0.1;
    minhz_ind = min(find(freq_ax>minhz));

    keyboard
    
    %% Detect Mains components
    if S.do_mains
        spec = abs_ft(1:num_ics,1:floor(length(abs_ft)/2));
        spec(:,1:minhz_ind) = 0;
        [~, ind] = max(spec');
        fmax = freq_ax(ind);
        keyboard
        kurt = kurtosis(tc,[],2);
        kurt = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis

        mains_ind=find(S.mains_frequency-1<fmax & fmax<S.mains_frequency+1 & kurt'<S.mains_kurt_thresh);
        
        if(length(mains_ind) > S.max_num_artefact_comps),
            warning([num2str(length(mains_ind)) ' Mains comps have been detected. Will only label the top ' num2str(S.max_num_artefact_comps) ' with the lowest kurtosis as artefacts.' ]);
            [~, new_inds] = sort(kurt(mains_ind));
            mains_ind = mains_ind(new_inds(1:S.max_num_artefact_comps));
        end
        
        if ~isempty(mains_ind)
            for j = 1:length(mains_ind)
                fprintf('Rejecting IC %d due to mains at %d Hz (Kurtosis=%.2f)\n',mains_ind(j),S.mains_frequency,kurt(mains_ind(j)));
                artefacts{mains_ind(j)} = sprintf('%s Mains',artefacts{mains_ind(j)});
            end
        else
            fprintf('No ICs associated with mains\n')
        end
    end

    % CHANNEL CORRELATIONS
    if(~isempty(S.artefact_channels) && length(S.artefact_channels)>0)
        artefact_chantype = unique(S.artefact_channels);
        for j = 1:length(artefact_chantype)
            reject = find(metrics.(artefact_chantype{j}).value > S.artefact_chans_corr_thresh(j));
            for k = 1:length(reject)
                fprintf('Rejecting IC %d due to %s (correlation = %.2f)\n',reject(k),artefact_chantype{j},metrics.(artefact_chantype{j}).value(reject(k)));
                artefacts{reject(k)} = sprintf('%s %s',artefacts{reject(k)},artefact_chantype{j});
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Artefact rejection using kurtosis.

    reject_kurt = [];
    if S.do_kurt > 0,
        [comps2reject,fig_handles_tmp,fig_names_tmp,fig_titles_tmp] = rank_and_plot(tc,D.ica.sm,abs_ft,freq_ax,D,'abs_kurtosis',modalities,samples_of_interest,S.kurtosis_wthresh,S.kurtosis_thresh,S.max_num_artefact_comps,S.do_plots);
        figs.handles = [figs.handles fig_handles_tmp];
        figs.names = [figs.names fig_names_tmp];
        figs.titles = [figs.titles fig_titles_tmp];
        if ~isempty(comps2reject)
            for j = 1:length(comps2reject)
                fprintf('Rejecting IC %d due to kurtosis\n',comps2reject(j));
                artefacts{comps2reject(j)} = sprintf('%s %s',artefacts{comps2reject(j)},'Kurtosis')
            end
        else
            fprintf('No components rejected due to kurtosis\n');
        end
    end

    % Compute the bad components
    bad_components = find(~cellfun(@isempty,artefacts));

    % Plot if required
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rank_and_plot
% Subfunction to rank tICs by some summary metric and then manually flag
% artefact components from topoplots, time course and spectrum.

function [comps2reject, fig_handles, fig_names, fig_titles] = rank_and_plot(tc,sm,abs_ft,freq_ax,D,metric,modalities,samples_of_interest,wthresh,thresh,max_num_artefact_comps,do_plots)

    fig_handles=[];
    fig_names={};
    fig_titles={};

    cc=1;

    switch lower(metric)
        case 'kurtosis'
            met=kurtosis(tc,[],2)-3;
            figlab='Kurtosis';
            direction='descend';
            thresh=thresh;
        case 'log_kurtosis'
            met=log(kurtosis(tc,[],2)-3);
            figlab='log(Kurtosis)'; 
            direction='descend';
            thresh=log(thresh);
        case 'abs_kurtosis'
            met = kurtosis(tc,[],2);
            met = abs(demean(boxcox1(met))); % make distribution more normal to better balance low/high kurtosis
            figlab='abs(Kurtosis)'; 
            direction='descend';
            thresh=abs(thresh);
        otherwise
    end

    [met_ord, comp_ind_ordered_by_met]=sort(met,direction);

    % use robust glm to identify extreme kurtosis

    [comps2reject_ordered,fig_handle,fig_name,fig_title]=find_outliers(met_ord,wthresh,thresh,do_plots,figlab,max_num_artefact_comps);

    comps2reject=comp_ind_ordered_by_met(comps2reject_ordered);

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [outlier_inds,fig_handle,fig_name,fig_title]=find_outliers(data,w_thresh,data_thresh,do_plots,fig_label,max_num_outliers)
    fig_handle = [];
    fig_name = [];
    fig_title = [];

    dm=ones(length(data),1);
    [~,stats] = robustfit(dm,data,'bisquare',4.685,'off');

    outlier_inds1=1:length(data);
    if(w_thresh>0),
        outlier_inds1=find(stats.w<w_thresh);
    end

    outlier_inds2=1:length(data);
    if(~isempty(data_thresh)),
        outlier_inds2=find(data>data_thresh);
    end

    outlier_inds=intersect(outlier_inds1,outlier_inds2);

    [~,ia]=sort(stats.w(outlier_inds));

    outlier_inds_sorted=outlier_inds(ia);

    num_comps=length(outlier_inds_sorted);

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
end


