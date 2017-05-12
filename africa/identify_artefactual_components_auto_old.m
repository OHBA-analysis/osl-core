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
% HL+MWW 2013

function [bad_components, metrics,fig_handles, fig_names, fig_titles] = identify_artefactual_components_auto(D,S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set-Up
arg = inputParser;

arg.addParameter('max_num_artefact_comps',10);
arg.addParameter('do_plots',1);
arg.addParameter('modality',[]);

% Mains defaults
arg.addParameter('do_mains',false); 
arg.addParameter('mains_frequency',50); 
arg.addParameter('mains_kurt_thresh',0.4); 

% Kurtosis defaults
arg.addParameter('do_kurt',false); 
arg.addParameter('kurtosis_thresh',20); 
arg.addParameter('kurtosis_wthresh',0); 

% Artefact chans
arg.addParameter('artefact_chans_corr_thresh',0.2);

arg.parse(S.ident_params);
S.ident_params = arg.Results;

if ~isempty(S.artefact_channels)
    S.ident_params.artefact_chans_corr_thresh = ones(size(S.artefact_channels))*S.ident_params.artefact_chans_corr_thresh;
end

D = D.montage('switch',0);

modalities = unique(D.chantype(find(strncmpi(S.modality,D.chantype,3))));  % changed by DM

fig_handles = [];
fig_names   = [];
fig_titles  = [];

sm = D.ica.sm;
tc = (D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))').';

%% Select only good data for classification
samples_of_interest = ~all(badsamples(D,':',':',':'));
samples_of_interest = reshape(samples_of_interest,1,D.nsamples*D.ntrials);

tc_nans = tc;
tc_nans(:,~samples_of_interest) = NaN;
tc = tc(:,samples_of_interest );

if strcmp(S.modality,'EEG')   % changed by DM
    chan_inds=setdiff(find(any([strcmp(D.chantype,'EEG')],1)),D.badchannels);
    map_inds(find(any([strcmp(D.chantype,'EEG')],1))) = 1:numel(find(any([strcmp(D.chantype,'EEG')],1)));
else
    chan_inds=setdiff(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)),D.badchannels);
    map_inds(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1))) = 1:numel(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)));
end
%%


num_ics = D.ica.params.num_ics;
abs_ft  = abs(fft(demean(tc(:,:),2),[],2));
freq_ax = 0:(D.fsample/2)/(floor(length(abs_ft)/2)-1):(D.fsample/2);

minhz = 0.1;
minhz_ind = min(find(freq_ax>minhz));

confirmed_artefact_names = [];
confirmed_artefact_inds  = [];
metrics_names = [];
metrics = [];
max_num_artefact_comps = S.ident_params.max_num_artefact_comps; % max number of components that will be allowed to be labelled as bad in each category

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detect Mains components

if S.ident_params.do_mains

    spec = abs_ft(1:num_ics,1:floor(length(abs_ft)/2));
    spec(:,1:minhz_ind) = 0;
    [~, ind] = max(spec');
    fmax = freq_ax(ind);
    
    kurt = kurtosis(tc,[],2);
    kurt = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis
        
    mains_ind=find(S.ident_params.mains_frequency-1<fmax & fmax<S.ident_params.mains_frequency+1 & kurt'<S.ident_params.mains_kurt_thresh);
    
    if(length(mains_ind) > max_num_artefact_comps),
        warning([num2str(length(mains_ind)) ' Mains comps have been detected. Will only label the top ' num2str(max_num_artefact_comps) ' with the lowest kurtosis as artefacts.' ]);
        [~, new_inds] = sort(kurt(mains_ind));
        mains_ind = mains_ind(new_inds(1:max_num_artefact_comps));
    end;
    
    %figure; plot(fmax);figure; plot(freq_ax,spec);
    
    msg = sprintf('\n%s%d%s\n%',['Mains (' num2str(S.ident_params.mains_frequency) ') interference split over '], numel(mains_ind), ' components.');
    fprintf(msg);

    confirmed_artefact_inds = mains_ind;
    metrics = kurt(mains_ind);
    
    for jj=1:length(mains_ind),
        confirmed_artefact_names{length(confirmed_artefact_names)+1} = 'Mains';
        metrics_names{length(metrics_names)+1} = 'Kurt';
                
    end;            
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detect artefact channels related components

if(~isempty(S.artefact_channels) && length(S.artefact_channels)>0)

    artefact_chan_ind_list         = [];
    artefact_chan_names            = [];
    artefact_chans_corr_thresh_new = [];
    
    %% find all chan indexes for passed in artefact channel names
    for ii=1:length(S.artefact_channels),
        
        inds = find(strcmp(D.chantype,S.artefact_channels{ii}));
        
        if isempty(inds),
            % allow, as a backup, specifying channel exact name
            inds = find(strcmpi(chanlabels(D), S.artefact_channels{ii})); 
            if isempty(inds),
                disp(['Artefact channel ' S.artefact_channels{ii} ' not found.']);
            else
                inds = [];
            end
        end
        
        artefact_chan_ind_list = [artefact_chan_ind_list inds];
        
        for jj=1:length(inds),
            artefact_chan_names{length(artefact_chan_names)+1} = S.artefact_channels{ii};
            artefact_chans_corr_thresh_new(length(artefact_chans_corr_thresh_new)+1) = S.ident_params.artefact_chans_corr_thresh(ii);
        end
    end
    
    %% loop through artefact channels detecting correlated IC time courses
    for ii=1:length(artefact_chan_ind_list),
           
        artefact_tc = D(artefact_chan_ind_list(ii),:,:);
        if D.ntrials > 1
            artefact_tc = reshape(artefact_tc,1,[]);
        end 
        
        if(strcmp(artefact_chan_names{ii},'ECG'))
            artefact_tc = bandpass(artefact_tc,[1 48],D.fsample);
        end
        
        artefact_tc=artefact_tc.^2;
        
        artefact_tc_corr=zeros(num_ics,1);
        for kk=1:num_ics
            sig=tc(kk,:);
            %psig=smooth(sig.^2);
            psig=(sig.^2)';
            artefact_tc_corr(kk)=abs(corr(psig,artefact_tc(samples_of_interest)'));
        end

        artefact_ic_inds = find(artefact_tc_corr>artefact_chans_corr_thresh_new(ii));

        msg = [artefact_chan_names{ii} ' artefacts split over ', num2str(numel(artefact_ic_inds)), ' components.'];   
        disp(msg);

        if(length(artefact_ic_inds)>max_num_artefact_comps),
            warning([num2str(length(artefact_ic_inds)) ' ' artefact_chan_names{ii} ' comps have been detected. Will only label the top ' num2str(max_num_artefact_comps) ' with the highest correlation with the artefact channel as artefacts.' ]);
            [sorted_artefact_tc_corr, new_inds]=sort(artefact_tc_corr(artefact_ic_inds),'descend');
            artefact_ic_inds=artefact_ic_inds(new_inds(1:max_num_artefact_comps));
        end;
        
        confirmed_artefact_inds=[confirmed_artefact_inds(:) ; artefact_ic_inds(:)];
        metrics=[metrics(:) ; artefact_tc_corr(artefact_ic_inds(:))];
        
        for jj=1:length(artefact_ic_inds),
            confirmed_artefact_names{length(confirmed_artefact_names)+1}=artefact_chan_names{ii};            
            metrics_names{length(metrics_names)+1}='Corr';
        end;   
    end
       
end

[confirmed_artefact_inds,inds]=unique(confirmed_artefact_inds(:));

for kk=1:length(inds),confirmed_artefact_names_new{kk}=confirmed_artefact_names{inds(kk)};end;
if length(inds)>0, confirmed_artefact_names=confirmed_artefact_names_new; end;

for kk=1:length(inds),metrics_names_new{kk}=metrics_names{inds(kk)};end;
if length(inds)>0, metrics_names=metrics_names_new; end;

for kk=1:length(inds),metrics_new(kk)=metrics(inds(kk));end;
if length(inds)>0, metrics=metrics_new; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do plots for mains and chan artefact related ICs

if S.ident_params.do_plots
    if ~isempty(confirmed_artefact_inds)
        for ii=1:length(confirmed_artefact_inds),            
            
            c=1;
            ncols=max(numel(modalities),2);
            
            fig_handles_tmp=sfigure; set(fig_handles_tmp,'Position',[1 1 1500 1000]);
            fig_handles=[fig_handles,fig_handles_tmp];
            fig_names{length(fig_names)+1}=[confirmed_artefact_names{ii} '_ic' num2str(confirmed_artefact_inds(ii))];               
            fig_titles{length(fig_titles)+1}=['AFRICA: ' confirmed_artefact_names{ii} ', ic' num2str(confirmed_artefact_inds(ii))];               

            for m =1:numel(modalities)
                subplot(2,ncols,c); 
                component_topoplot(D,sm(:,confirmed_artefact_inds(ii)),modalities{m},true);
                title({['Component ' num2str(confirmed_artefact_inds(ii)) ': ' confirmed_artefact_names{ii} ', ' metrics_names{ii} '=' num2str(metrics(ii))] modalities{m}}); axis tight;
                %title(['IC' num2str(confirmed_artefact_inds(ii)) ', ' confirmed_artefact_names{ii} ' Artefacts' modalities(m)]);%axis tight;
                c=c+1;
            end
            
            c=1;
            
            subplot(2,ncols,ncols+c);
            plot(freq_ax,abs_ft(confirmed_artefact_inds(ii),1:floor(length(abs_ft)/2)));title([confirmed_artefact_names{ii} ' Artefacts : Frequency Spectrum']);xlabel('Frequency (Hz)');ylabel('Spectrum');axis tight;
            c=c+1;

            subplot(2,ncols,ncols+c);
           
            if D.ntrials == 1
                t = D.time;
            else
                t = (1:length(samples_of_interest))./D.fsample;
            end
            
            plot(t,tc_nans(confirmed_artefact_inds(ii),:));title([confirmed_artefact_names{ii} ' Artefacts']);xlabel('Time (s)');axis tight;
            
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Artefact rejection using kurtosis.

reject_kurt = [];
if S.ident_params.do_kurt > 0,
    [comps2reject,fig_handles_tmp,fig_names_tmp,fig_titles_tmp] = rank_and_plot(tc_nans,sm,abs_ft,freq_ax,D,'abs_kurtosis',modalities,samples_of_interest,S.ident_params.kurtosis_wthresh,S.ident_params.kurtosis_thresh,max_num_artefact_comps);
    fig_handles=[fig_handles, fig_handles_tmp];
    for ii=1:length(fig_names_tmp)
        fig_names{length(fig_names)+1}=fig_names_tmp{ii};
    end
    for ii=1:length(fig_titles_tmp)
        fig_titles{length(fig_titles)+1}=fig_titles_tmp{ii};
    end
    reject_kurt = [reject_kurt comps2reject]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT BAD COMPONENTS

bad_components=unique([confirmed_artefact_inds(:); reject_kurt(:)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rank_and_plot
% Subfunction to rank tICs by some summary metric and then manually flag
% artefact components from topoplots, time course and spectrum.

function [comps2reject, fig_handles, fig_names, fig_titles] = rank_and_plot(tc,sm,abs_ft,freq_ax,D,metric,modalities,samples_of_interest,wthresh,thresh,max_num_artefact_comps)

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

[comps2reject_ordered,fig_handle,fig_name,fig_title]=find_outliers(met_ord,wthresh,thresh,1,figlab,max_num_artefact_comps);

comps2reject=comp_ind_ordered_by_met(comps2reject_ordered);

fig_handles(cc)=fig_handle;
fig_names{cc}=fig_name;
fig_titles{cc}=fig_title;

cc=cc+1;

% if strcmp(direction,'descend'),
%     num_comps2reject=max(comps2reject);
% else
%     num_comps2reject=min(comps2reject);
% end;
% 
% num_comps=length(comps2reject);
%         
% if(num_comps>max_num_artefact_comps),
%     warning(['AFRICA detected ' num2str(num_comps) ' ICs with ' figlab ' > ' num2str(thresh) '. Only rejecting the ' num2str(max_num_artefact_comps) ' worst components.']);
% 
%     num_comps=max_num_artefact_comps;
% end;
% 
% comps2reject=comps2reject(1:num_comps);

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

edit
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [outlier_inds,fig_handle,fig_name,fig_title]=find_outliers(data,w_thresh,data_thresh,do_plots,fig_label,max_num_outliers)

dm=ones(length(data),1);
[~,stats] = robustfit(dm,data,'bisquare',4.685,'off');

outlier_inds1=1:length(data);
if(w_thresh>0),
    outlier_inds1=find(stats.w<w_thresh);
end;

outlier_inds2=1:length(data);
if(~isempty(data_thresh)),
    outlier_inds2=find(data>data_thresh);
end;

outlier_inds=intersect(outlier_inds1,outlier_inds2);

[~,ia]=sort(stats.w(outlier_inds));

outlier_inds_sorted=outlier_inds(ia);

num_comps=length(outlier_inds_sorted);
         
if(num_comps>max_num_outliers),
     warning(['AFRICA detected ' num2str(num_comps) ' ' fig_label ' ICs. Only rejecting the ' num2str(max_num_outliers) ' worst components.']);
 
     num_comps=max_num_outliers;
end;
 
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


