%% Introduction to MEG Work Experience programme
% This example script is designed to be used as part of a guided 1-hour overview of MEG analysis aimed
% at audiences with no/limited experience in neuroimaging (e.g. public). It is intended
% to be followed together with a presenter from OHBA who will cover complementary 
% background content. The basic outline of the session is 
% 
% 1. Review of typical task design
% 2. Common analysis in sensor space
% 3. Overview of coregistration and beamforming
% 4. Spatial maps and activity timecourse in source space

%%
% 

%% Review of experimental task design
% Instructor to cover

%% Trial response, individual trials and averaging
%%
% First, we are going to load in the data files
abspath = @(s) fullfile(osldir,'example_data','intro_to_meg',s); % Get files relative to this path
oat = osl_load_oat(abspath('sensorspace_erf.oat'));
D = spm_eeg_load(abspath('meg_recording'));
response = @(sensor,trial) plot(D.time,mean(D(sensor,:,trial),3));
average_response = @(sensor) plot(D.time,mean(D(sensor,:,:),3));

%% 
% The variable |D| holds the MEG data. Have a look at what's in it
D

%%
% We can look at the response
response(10,5)

%%
% This means plot the signal from the 10th sensor, and the 5th trial.
% Note that there are three types of sensors
D.sensors('MEG').chantype(1:10)

%% 
% Where are they?
pos = D.sensors('MEG').chanpos;
figure('Color','k','InvertHardCopy','off')
scatter3(pos(:,1),pos(:,2),pos(:,3),120,'b','filled')
hold on
for j = 1:3:size(pos,1)
	text(pos(j,1),pos(j,2),pos(j,3),num2str(j),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'Color','w','FontWeight','bold');
end
axis equal
axis vis3d
axis off

%%
% So each numbers on the plot corresponds to the first of three sensors at
% that location
%
% *ACTIVITIES*
%
% * What does the response look like at other sensors?
% * How much variability do you think there is in the response?

%%
% We can also look at the response, averaged over trials
figure
average_response(10)

%%
% Have a look at the average response for some different sensors
%
% *ACTIVITIES*
%
% * Why do you think we normally average over trials?
% * What would be the effect of averaging over sensors as well? Do you think
%   that would be a good idea?

%% Stats - COPEs and t-stats_
% Instructor to cover

%% Computing sensor stats


% Now we are ready to look at some statistics. First, pick the one with the
% strongest response
first_level_results=oat.first_level.results_fnames{1};
first_level_results=oat_load_results(oat,first_level_results); 
first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts 
S2=oat.first_level.report; 
S2.stats=first_level_results; 
S2.modality=oat.first_level.report.modality_to_do; 
[vox_ind_max time_ind_max freq_ind_max stats max_stat] = oat_find_max_stats( S2 ); 

%%
% Then, make a plot with it
S2=[]; 
S2.stats=stats; 
S2.chanlabel=stats.chanlabels{vox_ind_max}; % can change to this setting to plot the time courses at the voxel with the max t-stat 
S2.first_level_cons_to_do=first_level_cons_to_do; % plots all of these contrasts 
[fig_handle fig_name fig_title] = oat_plot_vox_stats(S2); 

%%
% *ACTIVITIES*
%
% * If you change 'vox_ind_max' to a different number, how do the plots
%   change?

%% Spatial variation in sensor task response
% Now, we want a more systematic way to analyze the differences in space.
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGPLANAR'; % can also set this to 'MEGMAG'
S2.first_level_contrast=[3]; % view faces-motorbikes contrast
S2.view_cope=1; % set to 0 to see the t-stat
[cfg, dats, fig_handle]=oat_stats_multiplotER(S2);

%%
% This plot is interactive
% * You can select sensors by drawing a box on the plot
% * You can then click inside the box to generate a time series obtained by averaging over
% all of the plots.
% * Finally, you can draw a box to select times on the plot, which will then be displayed
% on a 'topoplot' showing how they change over space
%
% *ACTIVITIES*
%
% * Change |S2.first_level_contrast| to a different contrast
% * Change |S2.view_cope| to 0 to see the t-statistic instead

%% Source space analysis: Concepts
% Instructor to introduce

%% MRI structural scans
% Have a look at the whole head

runcmd('fsleyes %s',abspath('mri_scan_head.nii'))

%%
% <<osl_fsleyes_head.png>>

%%
% But usually, we're only interested in the brain

runcmd('fsleyes %s',abspath('mri_scan_brain.nii.gz'))

%%
% <<osl_fsleyes_brain.png>>

%% Coregistering MEG and MRI
% First problem is that we aren't too interested in the rest of the head - but we need to 
% know where the sensors are in this picture. That's what coregistration is
D = spm_eeg_load(abspath('rhino_example'))
rhino_display(D);

%% Beamforming
% Instructor to cover

%% Source-space spatial maps
% And then, look at spatial pattern of response

fsleyes(abspath('beamformer_erf.oat/session1_wholebrain_first_level_dir/cope1_2mm.nii.gz'))

%%
% <<osl_fsleyes_task.png>>

%%
% *ACTIVITIES*
%
% * 'cope1_2mm' corresponds to the first COPE. Try changing 'cope' to 'tstat' and also change the contrast
% e.g. 'tstat3_2mm'
% * 

%% Source space stats timecourses
% Now we are in source space, we can ask what the timecourses look like for locations in the brain rather than
% for sensors.
oat = osl_load_oat(abspath('beamformer_erf.oat'),'wholebrain_first_level')
mni_coord=[4,-82,-8]; % Visual Cortex Voxel
S2=[];
S2.vox_coord=mni_coord;
S2.stats=oat.first_level.results_fnames{1};
S2.oat=oat;
S2.first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);

%%
% This will bring up a the COPE and tstat estimates across time for a voxel in
% visual cortex. Note the prominant response around 100ms
% This corresponds to a point in Right Hemisphere Fusiform Cortex. Note that
% the 100ms response does not appear here, rather the later Face specific
% response is more dominant.
%
% *ACTIVITIES*
%
% * Try changing the voxel to 32,-64,-18
% * Use fsleyes to find where this voxel is


