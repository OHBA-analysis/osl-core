function [Dx fig_names fig_handles fig_titles]=osl_detect_badevent_v2(S)

% [Dx,fig_handles]=osl_detect_badevent_v2(S)
%
% MWW 2013

Dold=S.D;

fig_handles=[];
fig_names={};
fig_titles={};

index=0;

S2=[];
S2.D=S.D;
fname=fnamedat(Dold);
[pth nm]=fileparts(fname);
S2.outfile=[pth filesep 'S' nm '.mat'];
Dx=spm_eeg_copy(S2);

if ~isfield(S, 'plot_basename'),
    S.plot_basename='';
end;
if ~isfield(S, 'plot_basetitle'),
    S.plot_basetitle=S.plot_basename;
end;

for mm=1:length(S.modalities),

    modality=S.modalities{mm};
    
    disp(['Doing ' modality]); 
    
    % all good trials
    trials = Dx.indtrial(Dx.condlist,'good');

    chan_list=find(strcmp(chantype(Dx),modality));

    % all good chans
    chanind = setdiff(chan_list, Dx.badchannels);

    if ~isempty(Dx.badchannels)
        disp(['Already labelled bad chans=' mat2str(Dx.badchannels)]);
    end;

    max_iters=S.max_iter;
    iters=0;

    cont_chans=ones(length(S.outlier_measure_fns),1);
    cont_evs=ones(length(S.outlier_measure_fns),1);

    display_iters=0;
    if(display_iters)
    for ii=1:length(S.outlier_measure_fns),
        gcf{ii}=sfigure;
    end;
    end;

    while (sum(cont_evs)>0 || sum(cont_chans)>0),

        cont_chans=ones(length(S.outlier_measure_fns),1);
        cont_evs=ones(length(S.outlier_measure_fns),1);

        %% now find bad chans

        for ii=1:length(S.outlier_measure_fns),
            dat=Dx(chanind,:,trials);
            datchan=feval(S.outlier_measure_fns{ii},reshape(dat,size(dat,1),size(dat,2)*size(dat,3)),[],2);

            dm=ones(length(datchan),1);
            [b,stats] = robustfit(dm,datchan,'bisquare',4.685,'off');

            wthresh_chan=S.outlier_wthresh_chan(ii);
            [bad_chan]=find(stats.w<wthresh_chan);
            bad_chan_w=stats.w(bad_chan);
            [sorted_w iw]=sort(bad_chan_w);
            sorted_bad_chan=bad_chan(iw);
            
            if(display_iters)
                sfigure(gcf{ii});subplot(221);        
                [hs hsx]=hist(datchan,max(20,round(length(datchan)/10)));      
                bar(hsx,hs);
                title([S.outlier_measure_fns{ii} '(chan), ' modality]);
                a1=axis;
                subplot(223);
                plot(datchan,1:length(chanind),'*g');ho;
                plot(datchan(bad_chan),bad_chan,'or');
                hold off;
                a2=axis;
                axis([a1(1) a1(2) a2(3) a2(4) ]);
            end;

            if isfield(S,'max_bad_channels') && ~isempty(S.max_bad_channels),
                
                if (length(badchannels(Dx))+length(sorted_bad_chan)) > S.max_bad_channels,
                    warning(['More than S.max_bad_channels=' num2str(S.max_bad_channels) ' have been detected. But only marking the worst ' num2str(S.max_bad_channels) ' as bad']);
                end;
                    
                num_to_add=min(length(sorted_bad_chan),S.max_bad_channels-length(badchannels(Dx)));
                
                if num_to_add>0,
                    sorted_bad_chan=sorted_bad_chan(1:num_to_add);
                else
                    sorted_bad_chan=[];
                end;
            end;
            
            % set bad channels in Dx
            if(length(sorted_bad_chan)>0),                
                
                bad_channels=chanind(sorted_bad_chan);

                Dx = badchannels(Dx, bad_channels, ones(length(bad_channels),1));
                badchannels(Dx)
                cont_chans(ii)=1;

            else
                cont_chans(ii)=0;
            end;

            % correct chan inds for next bit
            chanind = setdiff(chan_list, Dx.badchannels);
        end;

        if(~S.just_chans)
            %% now find bad events
            for ii=1:length(S.outlier_measure_fns),
                dat=Dx(chanind,:,trials);
                datchan=feval(S.outlier_measure_fns{ii},reshape(dat,size(dat,1)*size(dat,2),size(dat,3)),[],1);

                dm=ones(length(trials),1);
                [b,stats] = robustfit(dm,datchan,'bisquare',4.685,'off');

                wthresh_ev=S.outlier_wthresh_ev(ii);
                bad_ev=find(stats.w<wthresh_ev);

                if(display_iters)
                    sfigure(gcf{ii});subplot(222);        
                    [hs hsx]=hist(datchan,max(20,round(length(trials)/10)));        
                    bar(hsx,hs);
                    title([S.outlier_measure_fns{ii} '(ev), ' modality]);
                    a1=axis;
                    subplot(224);
                    plot(datchan,1:length(trials),'*g');ho;
                    plot(datchan(bad_ev),bad_ev,'or');
                    hold off;
                    a2=axis;
                    axis([a1(1) a1(2) a2(3) a2(4) ]);
                end;

                if(length(bad_ev)>0),
                    bad_events=trials(bad_ev);
                    rejtmp=badtrials(Dx);
                    rej=zeros(1,Dx.ntrials);
                    rej(rejtmp)=1;
                    rej(bad_events)=1;
                    Dx = badtrials(Dx, 1:length(rej), rej);  

                    cont_evs(ii)=1;
                else
                    cont_evs(ii)=0;
                end;

                % correct trials for next bit
                trials = Dx.indtrial(Dx.condlist,'good');
            end;
        end;
        
        iters=iters+1;

        if(iters>max_iters)
            break;
        end;    

    end;

    if(S.do_plot),
        %% diagnostic plot
        pc=0.4;
        for ii=1:length(S.outlier_measure_fns),

            index=index+1;
            fig_handles(index)=sfigure;
            fig_names{index}=[S.plot_basename '_' S.outlier_measure_fns{ii} '_chans_' modality];            
            fig_titles{index}=[S.plot_basetitle S.outlier_measure_fns{ii} ', ' modality];
            
            %% chans

            % all good chans in orig
            chan_list=find(strcmp(chantype(Dold),modality));
            chanind = setdiff(chan_list, Dold.badchannels);

            chanindDx = setdiff(chan_list, Dx.badchannels);

            chanindDx_inchanind=zeros(length(chanindDx),1);
            for jj=1:length(chanindDx),
                tmp=find(chanind==chanindDx(jj));
                chanindDx_inchanind(jj)=tmp;            
            end;

            % all good trials in new
            trials = Dx.indtrial(Dx.condlist,'good');

            dat=Dold(chanind,:,trials);
            datDx=Dx(chanindDx,:,trials);

            datchan=feval(S.outlier_measure_fns{ii},reshape(dat,size(dat,1),size(dat,2)*size(dat,3)),[],2);        
            datchanDx=nan(size(datchan));
            datchanDx(chanindDx_inchanind)=feval(S.outlier_measure_fns{ii},reshape(datDx,size(datDx,1),size(datDx,2)*size(datDx,3)),[],2);

            snugplot(4,2,1,pc);
            [hs hsx]=hist(datchan,max(20,round(length(datchan)/10)));
            bar(hsx,hs);
            
            title([S.outlier_measure_fns{ii} '(chans), ' modality]);
            a1=axis;
            snugplot(4,2,3,pc);
            plot(datchan,1:length(chanind),'*r');ho;
            plot(datchanDx,1:length(chanind),'og');
            ylabel('Chan index');
            a2=axis;
            axis([a1(1) a1(2) a2(3) a2(4) ]);

            snugplot(4,2,5,pc);
            [hs hsx]=hist(datchanDx,max(20,round(length(datchanDx)/10)));
            bar(hsx,hs);
            a1=axis;
            snugplot(4,2,7,pc);
            plot(datchanDx,1:length(chanind),'*g');
            a2=axis;
            axis([a1(1) a1(2) a2(3) a2(4) ]);
            xlabel([S.outlier_measure_fns{ii} '(chans), ' modality]);
            ylabel('Chan index');
            
            if(~S.just_chans)
                %% evs

                % all good chans in new
                chan_list=find(strcmp(chantype(Dx),modality));
                chanind = setdiff(chan_list, Dx.badchannels);

                clear trials;
                % all good trials in passed in orig
                trials_old = Dold.indtrial(Dold.condlist,'good');
                trials_new = Dx.indtrial(Dx.condlist,'good');


                trials_new_intrials_old=zeros(length(trials_new),1);
                for jj=1:length(trials_new),
                    tmp=find(trials_old==trials_new(jj));
                    trials_new_intrials_old(jj)=tmp;            
                end;

                dat_old=Dold(chanind,:,trials_old);
                dat_new=Dx(chanind,:,trials_new);

                datchan_old=feval(S.outlier_measure_fns{ii},reshape(dat_old,size(dat_old,1)*size(dat_old,2),size(dat_old,3)),[],1);
                datchan_new=nan(size(datchan_old));
                datchan_new(trials_new_intrials_old)=feval(S.outlier_measure_fns{ii},reshape(dat_new,size(dat_new,1)*size(dat_new,2),size(dat_new,3)),[],1);

                snugplot(4,2,2,pc);
                [hs hsx]=hist(datchan_old,max(20,round(length(trials_old)/10)));
                bar(hsx,hs);
                title([S.outlier_measure_fns{ii} '(trials), ' modality]);
                a1=axis;
                snugplot(4,2,4,pc);
                plot(datchan_old,1:length(trials_old),'*r');ho;
                plot(datchan_new,1:length(trials_old),'og');
                a2=axis;
                axis([a1(1) a1(2) a2(3) a2(4) ]);
                ylabel('Trial index');
                
                snugplot(4,2,6,pc);
                [hs hsx]=hist(datchan_new,max(20,round(length(trials_new)/10)));
                bar(hsx,hs);
                a1=axis;
                snugplot(4,2,8,pc);
                plot(datchan_new,1:length(trials_old),'*g');ho;
                a2=axis;
                axis([a1(1) a1(2) a2(3) a2(4) ]);
                ylabel('Trial index');
                
                xlabel([S.outlier_measure_fns{ii} '(trials), ' modality]);
                
            end;
        end;     
    end;
end;

chanind = setdiff(chan_list, Dx.badchannels);

disp(['Bad chans=' mat2str(Dx.badchannels)]);

Dx.save;