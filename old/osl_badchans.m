%omt_badchans

function osl_badchans(fname)

D=spm_eeg_load(fname);%'/home/disk3/hluckhoo/data/spm_files/apoe_eo/spm8_apoe_eo001'

chans2plot=5;   % Plots five channels at a time. Do not adjust.
std_mult=3;     % Multiplier for STD plotting later. Can be adjusted. Isn't important.
ds=50;          % Downsampling factor. Must be integer. High values will speed things up but might overlook high frequency artefacts.
threshold=0.01; % Probalistic threshold for recommendation of good or bad. Not important at the moment but might be used later for automation.
fig_siz=[1200 600];
chan_deep=21;    % How many blocks of 5 channels are being plotted. 
%% Channel Grouping

mag_ind=find(strcmp('MEGMAG',D.chantype));
plan_ind=find(strcmp('MEGPLANAR',D.chantype));
mag_dat=D(mag_ind,:,:);
plan_dat=D(plan_ind,:,:);

var_mag=var(mag_dat,0,2); [varmag_ord mag_ord]=sort(var_mag,'descend');
var_plan=var(plan_dat,0,2);[varplan_ord plan_ord]=sort(var_plan,'descend');

%% Gamma Fitting: Channel Variances modelled as a gamma distribution.

bins_plan=0:0.005*max(varplan_ord):1.1*max(varplan_ord);bins_mag=0:0.0051*max(varmag_ord):1.1*max(varmag_ord);
hist_plan=hist(varplan_ord,bins_plan);hist_mag=hist(varmag_ord,bins_mag);

[phat_mag,pci] = gamfit(varmag_ord);
[phat_plan,pci] = gamfit(varplan_ord);
p_mag = pdf('Gamma',bins_mag,phat_mag(1),phat_mag(2));
p_plan = pdf('Gamma',bins_plan,phat_plan(1),phat_plan(2));

[phat_mag,pci] = gamfit(varmag_ord/std(varmag_ord));
[phat_plan,pci] = gamfit(varplan_ord/std(varplan_ord));
prob_var_mag=pdf('Gamma',varmag_ord/std(varmag_ord),phat_mag(1),phat_mag(2));
prob_var_plan = pdf('Gamma',varplan_ord/std(varplan_ord),phat_plan(1),phat_plan(2));

ymax_10_mag=(std_mult+2)*sqrt(mean(varmag_ord(1:10))); %ymax_90_mag=std_mult*sqrt(mean(varmag_ord(11:end)));
ymax_10_plan=(std_mult+2)*sqrt(mean(varplan_ord(1:20))); %ymax_90_plan=std_mult*sqrt(mean(varplan_ord(21:end)));

 [~, ind_mag]=min(abs(pdf('Gamma',bins_mag/std(varmag_ord),phat_mag(1),phat_mag(2))-0.05));
 ymax_90_mag=sqrt(bins_mag(ind_mag))*std_mult;
 [~, ind_plan]=min(abs(pdf('Gamma',bins_plan/std(varplan_ord),phat_plan(1),phat_plan(2))-0.05));
 ymax_90_plan=sqrt(bins_mag(ind_plan))*std_mult;
 
 
%% Visualisation
 
figure; 
subplot(1,2,1);bar(bins_mag, hist_mag); ho; plot(bins_mag,  max(hist_mag)*p_mag/max(p_mag),'r'); title('Magentometer Variance Distribution')
subplot(1,2,2);bar(bins_plan, hist_plan);ho; plot(bins_plan,  max(hist_plan)*p_plan/max(p_plan),'r'); title('Gradiometer Variance Distribution')

figure('Position',[100 100 fig_siz(1) fig_siz(2)])
disp('Looking at Magnetometers')
for j=1:chan_deep;
    for i=1:chans2plot
        if i+(j-1)*chans2plot<=102
            hax = axes('Position', [0.2, 1.05-0.2*i, 0.75, 0.1]);
            plot(downsample(D.time,ds), downsample(mag_dat(mag_ord(i+(j-1)*chans2plot),:),ds)); axis([min(D.time) max(D.time) -ymax_10_mag ymax_10_mag])
            hold on; plot(downsample(D.time,ds), ymax_90_mag*ones(1,length(downsample(D.time,ds))),'r--');plot(downsample(D.time,ds), -ymax_90_mag*ones(1,length(downsample(D.time,ds))),'r--');
            if prob_var_mag(i+(j-1)*chans2plot)<threshold; recommend = 'Recommend setting as bad.'; else recommend = 'Recommend leaving as good.'; end
            tit_h=title(['Magnetometer Channel ' D.chanlabels{mag_ind(mag_ord(i+(j-1)*chans2plot))} '; Variance is ' num2str(var_mag(mag_ord(i+(j-1)*chans2plot))) '. Probability equals ' num2str(prob_var_mag(i+(j-1)*chans2plot)) '. ' recommend]);
        end
    end
    ind=[1:5]+(j-1)*chans2plot; ind = ind(ind<=102);
    chans2consider=mag_ind(mag_ord(ind));
    selchan=input('Please flag bad channels for removal (e.g. [1 3 5]):     ');
    if ~isempty(selchan)
        chans_marked_as_bad=chans2consider(selchan);
        D=badchannels(D, D.meegchannels, ismember(D.meegchannels,[D.badchannels chans_marked_as_bad]));
        disp([{'Channels: '} D.chanlabels(find(ismember(D.meegchannels,[D.badchannels chans_marked_as_bad]))) {' are marked as bad'}])
    end
    clf;
end

disp('Looking at Gradiometers')
for j=1:2*chan_deep;
    for i=1:chans2plot
        if i+(j-1)*chans2plot<=204
        hax = axes('Position', [0.2, 1.05-0.2*i, 0.75, 0.1]);
        plot(downsample(D.time,ds), downsample(plan_dat(plan_ord(i+(j-1)*chans2plot),:),ds)); axis([min(D.time) max(D.time) -ymax_10_plan ymax_10_plan])
        hold on; plot(downsample(D.time,ds), ymax_90_plan*ones(1,length(downsample(D.time,ds))),'r--');plot(downsample(D.time,ds), -ymax_90_plan*ones(1,length(downsample(D.time,ds))),'r--');
        if prob_var_plan(i+(j-1)*chans2plot)<threshold; recommend = 'Recommend setting as bad.'; else recommend = 'Recommend leaving as good.'; end
        tit_h=title(['Gradiometer Channel ' D.chanlabels{plan_ind(plan_ord(i+(j-1)*chans2plot))} '; Variance is ' num2str(var_plan(plan_ord(i+(j-1)*chans2plot))) '. Probability equals ' num2str(prob_var_plan(i+(j-1)*chans2plot)) '. ' recommend]);
        end
    end
    ind=[1:5]+(j-1)*chans2plot; ind = ind(ind<=204);
    chans2consider=plan_ind(plan_ord(ind));
    selchan=input('Please flag bad channels for removal (e.g. [1 3 5]):     ');
    if ~isempty(selchan)
        chans_marked_as_bad=chans2consider(selchan);
        D=badchannels(D, D.meegchannels, ismember(D.meegchannels,[D.badchannels chans_marked_as_bad]));
        disp([{'Channels: '} D.chanlabels(find(ismember(D.meegchannels,[D.badchannels chans_marked_as_bad]))) {' are marked as bad'}])
        
    end
    clf;
end
D.badchannels;
save(D);
end
