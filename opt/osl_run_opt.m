function opt=osl_run_opt(opt)

% opt=osl_run_opt(opt)
%
% See osl_check_opt for settings
%
% See https://ohba-analysis.github.io/osl-docs/pages/docs/opt.html for
% documentation on OPT
%
% MWW 2013

opt = osl_check_opt(opt);

opt.results=[];

if(isempty(findstr(opt.dirname, '.opt')))
    opt.dirname=[opt.dirname, '.opt'];
end

mkdir(opt.dirname);

% set logfile up
opt.results.plotsdir=fullfile(opt.dirname, 'plots');
opt.results.logfile=fullfile(opt.results.plotsdir,['log-' date '.txt']);
mkdir(opt.results.plotsdir);

% delete any existing diary file with the same name
delete(opt.results.logfile);

opt.results.date=date;

% set diagnostic report up
opt_top_report=osl_report_setup(opt.results.plotsdir,['OPT report'],opt.results.logfile);

spm_files_basenames={};
spm_files_epoched_basenames={};

if(isempty(opt.sessions_to_do) || length(opt.sessions_to_do)==0)
    error('There are no sessions to analyse');
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['%%%%%%%%%%%%%%%%%%%%%%%  STARTING LOOP OVER SESSIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

for subi=1:length(opt.sessions_to_do)
    subnum=opt.sessions_to_do(subi);

    % single session results container:
    opt_results=[];

    % have individual log file for this subject
    opt_results.logfile=fullfile(opt.results.plotsdir,['log-' date '-session' num2str(subnum) '.txt']);
    
    % delete any existing diary file with the same name
    warning off;
    delete(opt_results.logfile)
    warning on;
    
    diary(opt_results.logfile);
    diary on;
    
    input_file=opt.spm_files{subnum};

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OPT ON SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  Input file: ' input_file '  %%%%%%%%%%%%%%%%%%%%%%%'])

    try

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Set session specific diagnostic report up
        report_dir=fullfile(opt.results.plotsdir,['session' num2str(subnum)]);
        opt_report=osl_report_setup(report_dir,['Session ' num2str(subnum) ' (Input file:' input_file  ')']);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load in passed in SPM MEEG object and setup filenames         
        D = spm_eeg_load(opt.spm_files{subnum});
        spm_files_basenames{subnum} = D.fname;
        
        Dnew = copy(D,[opt.dirname filesep D.fname]);

        % ensure there is no .mat extension:
        [p spmname e] = fileparts(spm_files_basenames{subnum});
        spm_files_basenames{subnum}=[spmname];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Downsample with SPM (particularly important if movement compensation used - but worth doing anyway)

        if(opt.downsample.do)
            S=[];
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.downsample, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            S.D=spm_file;
            S.fsample_new = opt.downsample.freq;
            D = spm_eeg_downsample (S);

            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                Dold.delete;
            end

            spm_files_basenames{subnum}=['d' spm_files_basenames{subnum}];

            close all
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% High Pass Filter
        % To get rid of low frequency fluctuations in the data

        if(opt.highpass.do)
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.highpass, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            S=[];
            S.D = spm_file;
            S.band='high';
            S.freq=opt.highpass.cutoff;
            S.dir='twopass';
            D = spm_eeg_filter(S);

            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                Dold.delete;
            end

            spm_files_basenames{subnum}=['f' spm_files_basenames{subnum}];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Mains Filter
        if(opt.mains.do)
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.mains, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
                       
            D = spm_eeg_load(spm_file); 
            
            if D.fsample>100           
                D = osl_filter(D,-1*(50+[-2 2])); % Remove 50Hz with notch filter
            else
                warning('Unable to do mains notch filtering at 50Hz, as D.fsample is not >100Hz');
            end
            
            if D.fsample>200            
                D = osl_filter(D,-1*(100+[-2 2])); % Remove 100Hz with notch filter
            else
                disp('Unable to do mains notch filtering at 100Hz, as D.fsample is not >200Hz');            
            end
            
            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                Dold.delete;
            end

            spm_files_basenames{subnum}=['f' spm_files_basenames{subnum}];
        end

        %%%%%%%%%%%%%%%%%%%
        %% DO REGISTRATION AND RUN FORWARD MODEL BASED ON STRUCTURAL SCANS
        if(opt.coreg.do)
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.coreg, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            S=[];
            spm_file=[opt.dirname filesep spm_files_basenames{subnum} '.mat'];
            S.D = [spm_file];
            S.mri=opt.coreg.mri{subnum};
            S.useheadshape=opt.coreg.useheadshape;
            S.forward_meg=opt.coreg.forward_meg;
            S.use_rhino=opt.coreg.use_rhino;

            if(isfield(opt.coreg,'fid_mnicoords')),
                S.fid.coords = opt.coreg.fid_mnicoords;
                S.fid.coordsys = 'MNI';
                % flirt -in /Users/woolrich/Desktop/GN170_anatomy_test.nii -ref /usr/local/fsl/data/standard/MNI152_T1_2mm -out /Users/woolrich/Desktop/anat_mne2;
            end
            S.fid.label = opt.coreg.fid_label;

            D=osl_headmodel(S);
            clc
            close all;

            % GENERATE REGISTRATION REPORTS

            spm_file=[opt.dirname filesep spm_files_basenames{subnum} '.mat'];
            Dcheck=spm_eeg_load(spm_file);
            
            %% spm displays from spm_eeg_inv_checkdatareg
            coregfig1 = sfigure;
            set(coregfig1,'Position',[300 300 1200 800]);
            set(coregfig1,'Color', [1,1,1]);
            subplot(1,3,1)
            warning off;spm_eeg_inv_checkdatareg_3Donly(Dcheck);warning on;
            view(-180,0)
            zoom(0.5);
            %title(['concatMefsession' num2str(counter) '_spm_meeg'])
            subplot(1,3,2)
            warning off;spm_eeg_inv_checkdatareg_3Donly(Dcheck);warning on;
            view(-270,0)
            zoom(0.5);
            subplot(1,3,3)
            warning off;spm_eeg_inv_checkdatareg_3Donly(Dcheck);warning on;
            view(130,18)
            zoom(0.5);
            %title(['Session ' num2str(subnum)]);

            opt_report=osl_report_set_figs(opt_report,'opt-coreg_spm_view',coregfig1);
            opt_report=osl_report_print_figs(opt_report);

            if opt.coreg.use_rhino
                %% now do rhino displays
                coregfig2 = sfigure;
                rhino_display(Dcheck,coregfig2);
                view(45,5)
                title(['Session ' num2str(subnum)]);
                opt_report=osl_report_set_figs(opt_report,'opt-coreg_rhino_view1',coregfig2);
                opt_report=osl_report_print_figs(opt_report);

                coregfig3 = sfigure;
                rhino_display(Dcheck,coregfig3);
                view(90,-10)
                %title(['Session ' num2str(subnum)]);
                opt_report=osl_report_set_figs(opt_report,'opt-coreg_rhino_view2',coregfig3);
                opt_report=osl_report_print_figs(opt_report);
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Perform AfRICA - ICA denoising

        if(opt.africa.todo.ica) || (opt.africa.todo.ident) || (opt.africa.todo.remove)

            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.africa, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            %% ADAPTER CODE FOR NEW CHANGES TO AFRICA
            S = struct;
            S.do_ica = opt.africa.todo.ica;
            S.do_ident = opt.africa.todo.ident;
            S.do_remove = opt.africa.todo.remove;
            S.precompute_topos = opt.africa.precompute_topos;
            S.used_maxfilter = opt.africa.used_maxfilter;
            S.mains_frequency = opt.africa.ident.mains_frequency;
            S.artefact_channels = opt.africa.ident.artefact_chans;
            S.auto_max_num_artefact_comps = opt.africa.ident.max_num_artefact_comps;
            S.auto_do_mains = opt.africa.ident.do_mains
            S.auto_mains_kurt_thresh = opt.africa.ident.mains_kurt_thresh
            S.auto_do_kurt = opt.africa.ident.do_kurt
            S.auto_kurtosis_thresh = opt.africa.ident.kurtosis_thresh
            S.auto_kurtosis_wthresh = opt.africa.ident.kurtosis_wthresh
            S.auto_artefact_chans_corr_thresh = opt.africa.ident.artefact_chans_corr_thresh

            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            D = spm_eeg_load(D);
            D=osl_africa(D,S);
            figs = report.ica(D);
            
            D.save()

            % If we removed artefacts, then we want to use the denoised data going forward
            % However, OAT does not work with online montages, so we need to clone 
            % the MEEG and store the AFRICA output on disk
            if S.do_remove
                [dir,nam,~] = fileparts(fullfile(D.path,D.fname));
                D2=clone(D.montage('switch',0),[dir '/A' nam '.dat'],size(D.montage('switch',0)));
                D2(:,:) = D(:,:);
                D2 = D2.montage('remove',1:D2.montage('getnumber'));
                D2.save();
                D = D2;
                spm_files_basenames{subnum}=['A' spm_files_basenames{subnum}];                    
            else
                warning('AFRICA has been run, but bad components have not been removed');
            end

            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) get(x,'tag'),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);

        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Mark bad segments
        if(opt.bad_segments.do)

            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.bad_segments, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])            
                            
            % copy file and add R prefix
            S=[];
            S.D=spm_file;
            S.outfile=prefix(spm_file,'B'); 
            D_continuous=spm_eeg_copy(S); 
            
            S = struct();
            S.modalities=opt.modalities;
            S.dummy_epoch_tsize    = opt.bad_segments.dummy_epoch_tsize;   
            S.measure_fns  = opt.bad_segments.outlier_measure_fns; 
            S.event_significance   = opt.bad_segments.event_significance;  
            S.channel_significance = opt.bad_segments.channel_significance;

            % first look for bad channels
            S.badchannels=true;
            S.badtimes=false;
            D_continuous = osl_detect_artefacts(D_continuous,S);
            D_continuous.save();

            S2=S;
            S2.plot_name_prefix='opt-bad_segments: ';
            figs = [report.bad_channels(D_continuous,opt.modalities, false, S2.plot_name_prefix)];
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) get(x,'tag'),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);
            
            % then look for bad segments
            S.badchannels=false;
            S.badtimes=true;
            D_continuous = osl_detect_artefacts(D_continuous,S);            
            D_continuous.save();

            S.plot_name_prefix='opt-bad_segments: ';
            figs = [report.bad_segments(D_continuous, S2)];
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) get(x,'tag'),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);

            ev = D_continuous.events;
            if isempty(ev)
                BadEpochs = [];
            else
                BadEpochs = ev(cellfun(@(x) ~isempty(strmatch('artefact',x)),{ev.type})); % Find only artefact events
            end
            opt_results.bad_segments.bad_segments=BadEpochs;
            
            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                Dold.delete;
            end

            spm_files_basenames{subnum}=['B' spm_files_basenames{subnum}];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %% Plot spectograms
        if(opt.spectograms.do)
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.spectograms, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            D_continuous=spm_eeg_load(spm_file);
            S.modalities=opt.modalities;

            for mm=1:length(S.modalities)
                fig_handles=sfigure;
                set(fig_handles,'Position',[1 1 1500 450]);

                S=[]; 
                S.D=D_continuous; 
                S.chantype=opt.modalities{mm};
                S.do_plot=false;
                S.cut_badsegments=true;
                [spect,F,T]=osl_plotspectrogram(S);

                %do plot
                imagesc(T,F,spect); 
                set(gca,'ydir','normal');
                plot4paper('time (s)','freq (hz)');

                opt_report=osl_report_set_figs(opt_report,['opt-spectogram_' opt.modalities{mm}],fig_handles,['opt-spectogram for ' opt.modalities{mm}]);
                opt_report=osl_report_print_figs(opt_report);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% DO EPOCHING (if epoch-based task data)

        if(opt.epoch.do)
            % DO preliminary epoching for the purpose of finding outliers
            % This is not the final epoching. Instead this sets up the epoch
            % definitions, and performs a temporary epoching for the purpose of doing
            % semi-automated outlier trial rejection (before running the fully
            % automated OAT).
            %
            % The epoch definitions and the continuous data will be kept and
            % passed into OAT. This is so that things like temporal filtering (which is
            % dones as part of OAT) can be done on the continuous data, before the data
            % is epoched inside OAT.
            %
            % Note that this will also remove those trials that overlap with the bad
            % epochs identified using oslview.

            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.epoch, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            %%%%
            % define the trials we want from the event information
            S2 = [];
            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];

            S2.D = spm_file;
            D_continuous=spm_eeg_load(S2.D);

            pretrig = opt.epoch.time_range(1)*1000; % epoch start in ms
            posttrig = opt.epoch.time_range(2)*1000; % epoch end in ms
            S2.timewin = [pretrig posttrig];

            S2.trialdef=opt.epoch.trialdef;
            S2.reviewtrials = 0;
            S2.save = 0;

            [epochinfo.trl, epochinfo.conditionlabels] = spm_eeg_definetrial(S2);

            %%%%
            % adjust timings to account for delay between trigger and visual display
            if(opt.epoch.timing_delay~=0)
                timing_delay=opt.epoch.timing_delay; % secs
                epochinfo_new=epochinfo;
                epochinfo_new.trl(:,1:2)=epochinfo.trl(:,1:2)+round(timing_delay/(1/D_continuous.fsample));
                epochinfo=epochinfo_new;
            end

            %%%%
            % do epoching
            S3 = epochinfo;
            S3.D = D_continuous;
            D = osl_epoch(S3);

            spm_files_epoched_basenames{subnum}=['e' spm_files_basenames{subnum}];

            opt_results.spm_files_epoched_basename=spm_files_epoched_basenames{subnum};
            epoched=true;
            
            % do plot of trial timings
            plot_name_prefix='opt-epoch: ';
            figs = report.trial_timings(D, [], plot_name_prefix);
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) sprintf('opt-epoch_%s',get(x,'tag')),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);            

        else
            % check if epoching already been done on spm_file

            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];

            D2=spm_eeg_load(spm_file);

            epoched= (D2.ntrials>1);

            if epoched
                spm_files_epoched_basenames{subnum}=spm_files_basenames{subnum};
                opt_results.spm_files_epoched_basename=spm_files_epoched_basenames{subnum};
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Detect bad trials and chans in epoched data        
        if(opt.outliers.do && epoched)

            disp(['%%%%%%%%%%%%%%%%%%%%%%%  opt.outliers, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            if ~opt.epoch.do           
                spm_files_epoched_basenames{subnum}=spm_files_basenames{subnum};
            end
            spm_file=[opt.dirname filesep spm_files_epoched_basenames{subnum}];

            % copy file for result of bad trial detection
            S2=[];
            S2.D=spm_file;
            %fname=fnamedat(spm_eeg_load(spm_file));
            [pth nm]=fileparts(spm_file);
            S2.outfile=[pth filesep 'R' nm '.mat'];
            D2=spm_eeg_copy(S2);
            
            % do bad trial/chan detection
            S = struct();
            S.modalities=opt.modalities;
            S.measure_fns  = opt.outliers.outlier_measure_fns;
            S.event_significance   = opt.outliers.event_significance;  
            S.channel_significance = opt.outliers.channel_significance;
            S.max_iter = 1;
            
            % first look for bad channels
            S.badchannels=true;
            S.badtimes=false;
            D2 = osl_detect_artefacts(D2,S);
            D2.save();
            
            plot_name_prefix='opt-outliers: ';            
            figs = [report.bad_channels(D2,opt.modalities, false, plot_name_prefix)];
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) sprintf('epoched_%s',get(x,'tag')),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);

            % then look for bad segments
            S.badchannels=false;
            S.badtimes=true;
            D2 = osl_detect_artefacts(D2,S);            
            D2.save();

            plot_name_prefix='opt-outliers: ';            
            figs = [report.bad_trials(D2, 'std', opt.modalities, plot_name_prefix)];
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) sprintf('epoched_%s',get(x,'tag')),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);

            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_epoched_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            
            % in case we had loaded in an SPM file from the start, check it
            % is not the exact same name as old SPM file?
            
            [p fname] = fileparts(opt.spm_files{subnum});
            if ((opt.cleanup_files == 1) || (opt.cleanup_files == 2)) && ~strcmp(Dold.fname, [fname '.mat'])
                Dold.delete;
            end

            spm_files_epoched_basenames{subnum}=['R' spm_files_epoched_basenames{subnum}];

            opt_results.spm_files_epoched_basename=spm_files_epoched_basenames{subnum};

        end

        opt_results.spm_files_basename=spm_files_basenames{subnum};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% write output filename
             
        if opt.epoch.do
            spm_file_result=[opt.dirname '/' opt_results.spm_files_epoched_basename '.mat'];   
        else
            spm_file_result=[opt.dirname '/' opt_results.spm_files_basename '.mat'];   
        end
        
        opt_report=osl_report_add_text(opt_report,['(Output file:' spm_file_result  ')'],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% write sub report and add it to top level report    
        opt_report=osl_report_write(opt_report, opt_top_report);
        opt_top_report=osl_report_add_sub_report(opt_top_report, opt_report);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save opt results
        opt_results.fname=['opt_result_session' num2str(subnum)];
        disp(['Saving opt results for: ' opt_results.fname]);
        opt_save_results(opt, opt_results);
        disp(['Save complete']);
        opt.results.fnames{subnum}=opt_results.fname;

    catch ME

         % set output SPM file as empty for this subject:
        spm_files_basenames{subnum}=[];
        if(opt.epoch.do)
            spm_files_epoched_basenames{subnum}=[];
        end

        disp(['Subject ' num2str(subnum) ' failed']);
        ME.getReport

    end
    
    % build logfile
    disp(['Building logfile']);
    diary(opt.results.logfile);
    text = fileread(opt_results.logfile);
    disp(text);

end

diary(opt.results.logfile);

%%%%%%%%%%%%%%%%%%%%
%% gather results over all sessions
disp(['Gathering results']);
opt=opt_gather_results(opt);

%%%%%%%%%%%%%%%%%%%%
%% diagnostic plots over all sessions
[opt, opt_top_report]=opt_report_summary_plots(opt, opt_top_report);

%%%%%%%%%%%%%%%%%%%
%% generate web report
disp(['Generating opt report']);
opt.results.report=osl_report_write(opt_top_report);

%%%%%%%%%%%%%%%%%%%
%% output res
opt.osl2_version=osl_version;

opt.fname=[opt.dirname filesep 'opt'];

disp(['Saving opt']);
save(opt.fname, 'opt');

disp(['To view OPT report, point your browser to <a href="' opt.results.report.html_fname '">' opt.results.report.html_fname '</a>']);

diary off;
