function opt=osl_run_opt(opt)

% opt=osl_run_opt(opt)
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
opt_report=osl_report_setup(opt.results.plotsdir,['OPT report'],opt.results.logfile);

spm_files_basenames={};
spm_files_epoched_basenames={};

if(isempty(opt.sessions_to_do) || length(opt.sessions_to_do)==0)
    error('There are no sessions to analyse');
end;

if opt.maxfilter.do && ~strcmp(opt.input_file_type,'raw_fif_files'),
    error('Input files must be raw file files to do Maxfiltering. Consider changing opt.maxfilter.do to 0.');
end;
if ~opt.maxfilter.do && strcmp(opt.input_file_type,'raw_fif_files'),
    error('Input files can not be raw fif files if NOT Maxfiltering. Consider changing opt.maxfilter.do to 1.');
end;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['%%%%%%%%%%%%%%%%%%%%%%%  STARTING LOOP OVER SESSIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

for subi=1:length(opt.sessions_to_do),
    subnum=opt.sessions_to_do(subi);

    % single session results container:
    opt_results=[];

    % have individual log file for this subject
    opt_results.logfile=fullfile(opt.results.plotsdir,['log-' date '-session' num2str(subnum) '.txt']);
    % delete any existing diary file with the same name
    delete(opt_results.logfile)

    diary(opt_results.logfile);
    diary on;
    
    switch opt.input_file_type
    case 'raw_fif_files'
        input_file=opt.raw_fif_files{subnum};
    case 'input_files'
        input_file=opt.input_files{subnum};
    case 'spm_files'
        input_file=opt.spm_files{subnum}
    otherwise
        error('Invalid input_file type');
    end;

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OPT ON SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  Input file: ' input_file '  %%%%%%%%%%%%%%%%%%%%%%%'])

    try

        % set session specific diagnostic report up
        report_dir=fullfile(opt.results.plotsdir,['session' num2str(subnum)]);
        opt_report=osl_report_setup(report_dir,['Session ' num2str(subnum) ' (Input file:' input_file  ')']);

        if(opt.maxfilter.do)

            opt_results.maxfilter=[];

            %if(~opt.maxfilter.do_sss_maxfilter ||  (opt.maxfilter.do_remove_badchans_pre_sss && opt.maxfilter.do_sss_maxfilter))

            %%%%%%%%%%%%%%%%%%
            %% Create 1st round of non-SSSed data if doing double Maxfilter
            % OR create final non-SSSed Maxfiltered data if not doing SSS
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  MAXFILT WITH NO SSS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            Smf=[];
            Smf.maxfilt_dir=opt.maxfilter.maxfilt_dir;
            [p fifname e] = fileparts(opt.raw_fif_files{subnum});
            Smf.fif=[p filesep fifname];

            if opt.maxfilter.ctc_file ~= 0
                Smf.ctc_file = ctc_file;
            end
            if opt.maxfilter.cal_file ~= 0
                Smf.cal_file = cal_file;
            end

            % output fif file from maxfilter call will be placed in the opt
            % directory with a name based on opt.convert.spm_files_basenames{subnum}
            fifname_out=['nosss_fif_' opt.convert.spm_files_basenames{subnum}];
            Smf.fif_out=[opt.dirname filesep fifname_out];
            Smf.logfile=1;

            if ~(exist([Smf.fif '.fif'],'file')==2),
                error(['opt.raw_fif_files{' num2str(subnum) '} file ' [Smf.fif '.fif'] ' does not exist']);
            end;

            Smf.movement_compensation=0;

            if(~Smf.movement_compensation)
                disp('Assuming data was sampled at 1000 Hz');
                sample_freq=1000;
                if(opt.downsample.freq<=sample_freq/4 && opt.downsample.do)
                    Smf.downsample_factor=4;
                elseif(opt.downsample.freq<=sample_freq/2 && opt.downsample.do)
                    Smf.downsample_factor=2;
                else
                    Smf.downsample_factor=1;
                end;
            end;

            Smf.bad_epochs=opt.maxfilter.bad_epochs{subnum};

            Smf.nosss=1;
            if opt.maxfilter.remote_port ~= 0
                fif_sss=osl_call_maxfilter_remote(Smf,opt.maxfilter.remote_port);
            else
                fif_sss=osl_call_maxfilter(Smf);
            end

            type([fif_sss '_log.txt']);

            if(~exist([fif_sss '.fif'],'file'))
                warning([fif_sss '.fif does not exist. Can not continue with this subject.']);
                break;
            end;

            %end;

            %%%%%%%%%%%%%%%%%%%
            %% Convert to SPM

            S2=[];
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  CONVERT, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            spm_file=fullfile(opt.dirname,['nosss' opt.convert.spm_files_basenames{subnum}]);
            S2.spm_file=spm_file;
            [p fifname e] = fileparts(Smf.fif_out);
            S2.fif_file=[p filesep fifname '.fif'];
            if(isfield(opt.convert,'trigger_channel_mask'))
                S2.trigger_channel_mask=opt.convert.trigger_channel_mask;
            end;
            [D] = osl_import(S2);

            % delete fif file in opt dir that is no longer needed
            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                delete(S2.fif_file);
            end;

            close all

            %%%%%%%%%%%%%%%%%%%
            %%
            if(opt.maxfilter.do_sss)

                if(opt.maxfilter.do_remove_badchans_pre_sss),

                    %%%%%%%%%%%%%%%%%%%
                    %% Do DOUBLE Maxfilter procedure:
                    %% First detects and remove bad channels on non-SSS data,
                    %% then runs SSS
                    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RM BAD CHANS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
                    spm_file=fullfile(opt.dirname ,['nosss' opt.convert.spm_files_basenames{subnum}]);
                    D=spm_eeg_load(spm_file);

                    for ii=1:length(opt.modalities),
                        printprefix_mod=[opt.modalities{ii} '_preSSS'];

                        %%%%%%%%%%%%%%%%%%%%%
                        % use outlier detection approach
                        if(1),

                            S = [];
                            S.D = D;
                            S.modalities{1}=opt.modalities{ii};
                            S.do_plot=0;
                            S.max_iter=1;
                            %S.outlier_measure_fns={'min','std'};
                            S.outlier_measure_fns={'std'};
                            S.outlier_wthresh_chan=[0.01];
                            S.just_chans=1;
                            S.plot_basename='double_maxfilter';
                            S.plot_basetitle='MAXFILTER: ';
                            S.max_bad_channels=opt.maxfilter.max_badchans_pre_sss;
                            [D2]=osl_detect_badevent(S);

                            list_bch=badchannels(D2);

                            if(~isempty(list_bch))
                                bad_list=badchannels(D);
                                vect_bch=zeros(1,size(D,1));
                                vect_bch(bad_list)=1;
                                vect_bch(list_bch)=1;

                                D = D.badchannels(cat(2,bad_list,list_bch),1);
                                D.save;
                            end;

                            % delete obsolete spm file
                            spm_file_old=[opt.dirname filesep 'Snosss' opt.convert.spm_files_basenames{subnum}];
                            Dold=spm_eeg_load(spm_file_old);

                            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                                Dold.delete;
                            end;

                            %%%%%%%%%%%%%
                            % diagnostic plots
                            S=[];
                            S.D=D;
                            S.print_plots=1;
                            S.plotsdir=opt.results.plotsdir;
                            S.modality=opt.modalities{ii};
                            S.printprefix=[printprefix_mod '_outlier'];
                            S.plot_basename='double_maxfilter';
                            S.plot_basetitle='MAXFILTER: ';
                            [res fig_names fig_handles fig_titles]=osl_check_bad_chans(S); % TODO - remove this call 

                            opt_report=osl_report_set_figs(opt_report,fig_names,fig_handles,fig_titles);
                            opt_report=osl_report_print_figs(opt_report);

                            %%%%%%%%%%%%%%%
                        end;
                        %%%%%%%%%%%%%%%%%%%%%

                        D.save;

                    end;

                    opt_results.maxfilter.badchans_pre_sss=badchannels(D);

                end;

                %%%%%%%%%%%%%%%%%%%
                %% Re-Maxfilter data using bad channel information

                disp(['%%%%%%%%%%%%%%%%%%%%%%%  MAXFILT WITH SSS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
                Smf=[];
                [p fifname e] = fileparts(opt.raw_fif_files{subnum});
                Smf.fif=[p filesep fifname];

                if opt.maxfilter.ctc_file ~= 0
                    Smf.ctc_file = ctc_file;
                end
                if opt.maxfilter.cal_file ~= 0
                    Smf.cal_file = cal_file;
                end

                % output fif file from maxfilter call will be placed in the opt
                % directory with a name based on opt.convert.spm_files_basenames{subnum}
                fifname_out=['sss_fif_' opt.convert.spm_files_basenames{subnum}];
                Smf.fif_out=[opt.dirname filesep fifname_out];

                Smf.logfile=1;
                Smf.nosss=0;

                spm_file=[opt.dirname filesep 'nosss' opt.convert.spm_files_basenames{subnum}]
                Smf.spmfile=spm_file;

                if ~isempty(opt.maxfilter.trans_ref_file) && ~strcmp(opt.maxfilter.trans_ref_file,''),
                    [pth fifname ext] = fileparts(opt.maxfilter.trans_ref_file);
                    Smf.trans_ref_file=[pth filesep fifname];
                end;

                Smf.movement_compensation=opt.maxfilter.movement_compensation; %can only run with sss AND downsampling has to be switched off

                if(opt.downsample.do && Smf.movement_compensation)
                    disp('Can not do any downsampling during Maxfilter call when Smf.movement_compensation is on.');
                end;

                Smf.st=[];
                Smf.st.do=opt.maxfilter.temporal_extension;

                if(~Smf.movement_compensation)
                    disp('Assuming data was sampled at 1000 Hz');
                    sample_freq=1000;
                    if(opt.downsample.freq<=sample_freq/4 && opt.downsample.do)
                        Smf.downsample_factor=4;
                    elseif(opt.downsample.freq<=sample_freq/2 && opt.downsample.do)
                        Smf.downsample_factor=2;
                    else
                        Smf.downsample_factor=1;
                    end;
                end;

                Smf.bad_epochs=opt.maxfilter.bad_epochs{subnum};

                if opt.maxfilter.remote_port ~= 0
                    fif_sss=osl_call_maxfilter_remote(Smf,opt.maxfilter.remote_port);
                else
                    fif_sss=osl_call_maxfilter(Smf);
                end
                type([fif_sss '_log.txt']);

                opt_results.maxfilter.sss_autobad_off(subnum)=0;

                % check if Maxfilter has worked
                spm_files_basenames{subnum}=['sss' opt.convert.spm_files_basenames{subnum}];
                [maxfilter_failed D opt_report]=opt_maxfilter_check_output(opt, fif_sss, spm_files_basenames{subnum}, opt_report);

                if maxfilter_failed,
                    warning('Maxfilter has FAILED. Will retry with no bad channels set and autobad off!!! Results may be suspect.');

                    opt_results.maxfilter.sss_autobad_off(subnum)=1;
                    % set all chans good
                    D=spm_eeg_load(spm_file);
                    chans=D.indchantype('MEEG');
                    D = D.badchannels(1:size(chans),0);
                    D.save;

                    Smf.autobad_off=1;
                    if opt.maxfilter.remote_port ~= 0
                        fif_sss=osl_call_maxfilter_remote(Smf,opt.maxfilter.remote_port);
                    else
                        fif_sss=osl_call_maxfilter(Smf);
                    end
                    type([fif_sss '_log.txt']);

                    % check if Maxfilter has worked
                    [maxfilter_failed D opt_report]=opt_maxfilter_check_output(opt, fif_sss, spm_files_basenames{subnum}, opt_report);

                    if maxfilter_failed,
                        error('Maxfilter has still FAILED. Can not continue with this subject.');
                    end;
                end;

                %%%
                % delete obsolete nosss spm file with pre sss bad channels
                spm_file_old=[opt.dirname filesep 'nosss' opt.convert.spm_files_basenames{subnum}];
                Dold=spm_eeg_load(spm_file_old);
                if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                    Dold.delete;
                end;
                %%%

            end;

        else,

            if(~isempty(opt.input_files)),
                % Skipping Maxfilter calls - assume these have already been run,
                % or we are not working with Neuromag data,
                % and so go straight for convert call using opt.input_files

                %%%%%%%%%%%%%%%%%%%
                %% Convert to SPM
                disp(['%%%%%%%%%%%%%%%%%%%%%%%  CONVERT INPUT DATA, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

                S2=[];

                [p spmname e] = fileparts(opt.convert.spm_files_basenames{subnum});
                spm_files_basenames{subnum}=[spmname '.mat'];
                S2.outfile=[opt.dirname filesep spm_files_basenames{subnum}];

                [p fifname e] = fileparts(opt.input_files{subnum});

                if isempty(e)
                    switch opt.datatype,
                        case 'neuromag'
                            e='.fif';
                        case 'ctf'
                            e='.ds';
                        otherwise
                            error('Unknown opt.datatype');
                    end;
                end;

                if(isfield(opt.convert,'trigger_channel_mask'))
                    S2.trigger_channel_mask=opt.convert.trigger_channel_mask;
                end;
                D = osl_import([p filesep fifname e],S2);
                fig_handles = report.events(D);
               
                if ~isempty(fig_handles)
                    fig_titles={}; fig_titles{1}='CONVERT: histogram of trigger codes';
                    opt_report=osl_report_set_figs(opt_report,{get(fig_handles,'tag')},fig_handles,fig_titles);
                    opt_report=osl_report_print_figs(opt_report);
                end

                close all
            else

                % use passed in SPM MEEG objects           
                D = spm_eeg_load(opt.spm_files{subnum});
                spm_files_basenames{subnum} = D.fname;

                Dnew = copy(D,[opt.dirname filesep D.fname]);

            end;

        end;

        % ensure there is no .mat extension:
        [p spmname e] = fileparts(spm_files_basenames{subnum});
        spm_files_basenames{subnum}=[spmname];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Use any bad epochs passed in
        if ~isempty(opt.convert.bad_epochs{subnum}),
            BadEpochs=opt.convert.bad_epochs{subnum};

            BadEpochs2=cell(size(BadEpochs,1),1);
            for ee=1:size(BadEpochs,1),
                BadEpochs2{ee}=BadEpochs(ee,:);
                if BadEpochs2{ee}(1)==-1, BadEpochs2{ee}(1)=D.time(1); end;
                if BadEpochs2{ee}(2)==-1, BadEpochs2{ee}(2)=D.time(end); end;
            end;

            D=set_bad(D,BadEpochs2);
            D.save;

            S2=[];
            S2.outlier_measure_fns={'std'};
            S2.modalities=opt.modalities;
            S2.BadEpochs=BadEpochs2;
            [fig_handles fig_names fig_titles] = plot_bad(D, S2);
            opt_report=osl_report_set_figs(opt_report, fig_names, fig_handles, fig_titles);
            opt_report=osl_report_print_figs(opt_report);
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Downsample with SPM (particularly important if movement compensation used - but worth doing anyway)

        if(opt.downsample.do),
            S=[];
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  DOWNSAMP, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            S.D=spm_file;
            S.fsample_new = opt.downsample.freq;
            D = spm_eeg_downsample (S);

            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            if(opt.cleanup_files == 1) || (opt.cleanup_files == 2)
                Dold.delete;
            end;

            spm_files_basenames{subnum}=['d' spm_files_basenames{subnum}];

            close all
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% High Pass Filter
        % To get rid of low frequency fluctuations in the data

        if(opt.highpass.do)
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  HP FILT, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
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
            end;

            spm_files_basenames{subnum}=['f' spm_files_basenames{subnum}];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Mains Filter
        if(opt.mains.do)
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  MAINS FILT, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

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
            end;

            spm_files_basenames{subnum}=['f' spm_files_basenames{subnum}];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Perform AfRICA - ICA denoising

        if(opt.africa.todo.ica) || (opt.africa.todo.ident) || (opt.africa.todo.remove)

            disp(['%%%%%%%%%%%%%%%%%%%%%%%  AFRICA, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

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
            
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  MARK  SEGMENTS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            D_continuous=spm_eeg_load(spm_file);
            S = struct();
            S.modalities=opt.modalities;
            S.dummy_epoch_tsize    = opt.bad_segments.dummy_epoch_tsize   
            S.measure_fns  = opt.bad_segments.outlier_measure_fns 
            S.event_significance   = opt.bad_segments.event_significance  
            S.channel_significance = opt.bad_segments.channel_significance
            D_continuous = osl_detect_artefacts(D_continuous,S);
            D_continuous.save();

            figs = [report.bad_channels(D_continuous) report.bad_segments(D_continuous)];
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) get(x,'tag'),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);

            ev = D_continuous.events;
            if isempty(ev)
                BadEpochs = [];
            else
                BadEpochs = ev(cellfun(@(x) ~isempty(strmatch('artefact',x)),{ev.type})); % Find only artefact events
            end
            opt_results.bad_segments.bad_segments=BadEpochs;

        else
            % plot data
            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];
            D_continuous=spm_eeg_load(spm_file);
            S.outlier_measure_fns={'std', 'min'};
            S.modalities=opt.modalities;
            for ss=1:length(S.outlier_measure_fns),
                fig_handles=sfigure;
                set(fig_handles,'Position',[1 1 1300 450]);
                for mm=1:length(S.modalities),
                    chan_inds = find(strcmp(D_continuous.chantype, S.modalities{mm}));
                    subplot(2,1,mm);
                    tmpdat=squeeze(D_continuous(chan_inds,:,1));
                    plot(D_continuous.time, feval(S.outlier_measure_fns{ss}, tmpdat));
                    a=axis;
                    ylims=a(3:4);

                    plot4paper('time(s)',[S.outlier_measure_fns{ss} '(' S.modalities{mm} ')']);
                end;
                opt_report=osl_report_set_figs(opt_report,[S.outlier_measure_fns{ss} '_chans_with_bad_epochs_tc'],fig_handles,['DATA: ' S.outlier_measure_fns{ss} ' over chans']);
                opt_report=osl_report_print_figs(opt_report);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %% plot spectograms
        disp(['%%%%%%%%%%%%%%%%%%%%%%%  PLOT SPECTOGRAMS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

        D_continuous=spm_eeg_load(spm_file);
        S.modalities=opt.modalities;

        for mm=1:length(S.modalities),
            fig_handles=sfigure;
            set(fig_handles,'Position',[1 1 1300 450]);

            S=[]; 
            S.D=D_continuous; 
            S.chantype=opt.modalities{mm};
            S.do_plot=false;
            S.cut_badsegments=true;
            [spect,F,T]=osl_plotspectrogram(S);

            %do plot
            imagesc(T,F,spect); colorbar
            set(gca,'ydir','normal');
            plot4paper('time (s)','freq (hz)');

            opt_report=osl_report_set_figs(opt_report,['spectogram_' opt.modalities{mm}],fig_handles,['DATA: spectogram for ' opt.modalities{mm}]);
            opt_report=osl_report_print_figs(opt_report);
        end

        %%%%%%%%%%%%%%%%%%%
        %% DO REGISTRATION AND RUN FORWARD MODEL BASED ON STRUCTURAL SCANS
        % needs to be done before any montaging (e.g.in AFRICA) to ensure that
        % Neuromag gradiometer baseline correction is done correctly.
        if(opt.coreg.do),
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  COREG, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
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
            end;
            S.fid.label = opt.coreg.fid_label;

            D=osl_headmodel(S);
            clc
            close all;

            % CHECK REGISTRATION

            % mnifid = ft_transform_headshape(D.inv{1}.datareg.toMNI,D.inv{1}.datareg.fid_mri );mnifid.fid.pnt
            % opt.fid_mnicoords.nasion =[  1.3968   81.9389  -44.9899];opt.fid_mnicoords.lpa =[-83.3246  -20.2372  -68.1528];opt.fid_mnicoords.rpa = [83.9906  -19.5985  -65.6612];

            %%
            spm_file=[opt.dirname filesep spm_files_basenames{subnum} '.mat'];
            Dcheck=spm_eeg_load(spm_file);
            %spm_eeg_inv_checkdatareg(D);

            %% spm displays
            coregfig1 = sfigure;
            set(coregfig1,'Position',[300 300 1200 800]);
            set(coregfig1,'Color', [1,1,1]);
            subplot(1,3,1)
            warning off;spm_eeg_inv_checkdatareg_3Donly(Dcheck);warning on;
            view(-180,0)
            %title(['concatMefsession' num2str(counter) '_spm_meeg'])
            subplot(1,3,2)
            warning off;spm_eeg_inv_checkdatareg_3Donly(Dcheck);warning on;
            view(-270,0)
            subplot(1,3,3)
            warning off;spm_eeg_inv_checkdatareg_3Donly(Dcheck);warning on;
            view(130,18)
            title(['Session ' num2str(subnum)]);
            %%
            opt_report=osl_report_set_figs(opt_report,'Coregistration_spm_view',coregfig1);
            opt_report=osl_report_print_figs(opt_report);

            if opt.coreg.use_rhino
                %% now do rhino displays
                coregfig2 = sfigure;
                rhino_display(Dcheck,coregfig2);
                view(45,5)
                title(['Session ' num2str(subnum)]);
                opt_report=osl_report_set_figs(opt_report,'Coregistration_rhino_view1',coregfig2);
                opt_report=osl_report_print_figs(opt_report);

                coregfig3 = sfigure;
                rhino_display(Dcheck,coregfig3);
                view(90,-10)
                title(['Session ' num2str(subnum)]);
                opt_report=osl_report_set_figs(opt_report,'Coregistration_rhino_view2',coregfig3);
                opt_report=osl_report_print_figs(opt_report);
            end
                
            %%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% DO EPOCHING (if epoch-based task data)

        if(opt.epoch.do),
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

            disp(['%%%%%%%%%%%%%%%%%%%%%%%  EPOCH, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

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

            [epochinfo.trl, epochinfo.conditionlabels, S3] = spm_eeg_definetrial(S2);

            %%%%
            % adjust timings to account for delay between trigger and visual display
            if(opt.epoch.timing_delay~=0)
                timing_delay=opt.epoch.timing_delay; % secs
                epochinfo_new=epochinfo;
                epochinfo_new.trl(:,1:2)=epochinfo.trl(:,1:2)+round(timing_delay/(1/D_continuous.fsample));
                epochinfo=epochinfo_new;
            end;

            %%%%
            % do epoching
            S3=[];
            S3 = epochinfo;
            S3.D = D_continuous;
            D = osl_epoch(S3);

            spm_files_epoched_basenames{subnum}=['e' spm_files_basenames{subnum}];

            opt_results.spm_files_epoched_basename=spm_files_epoched_basenames{subnum};
            epoched=true;
        else
            % check if epoching already been done on spm_file

            spm_file=[opt.dirname filesep spm_files_basenames{subnum}];

            D2=spm_eeg_load(spm_file);

            epoched= (D2.ntrials>1) ;

            if epoched
                spm_files_epoched_basenames{subnum}=spm_files_basenames{subnum};
                opt_results.spm_files_epoched_basename=spm_files_epoched_basenames{subnum};
            end
        end

        if(opt.outliers.do && epoched),
            %%%%%%%%%%%%%%%%%%%%
            %% Detect bad events and chans

            disp(['%%%%%%%%%%%%%%%%%%%%%%%  BAD CHAN/EVENTS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

            if ~opt.epoch.do           
                spm_files_epoched_basenames{subnum}=spm_files_basenames{subnum};
            end;
            spm_file=[opt.dirname filesep spm_files_epoched_basenames{subnum}];

            S = struct();
            S.modalities=opt.modalities;
            S.measure_fns  = opt.outliers.outlier_measure_fns 
            S.event_significance   = opt.outliers.event_significance  
            S.channel_significance = opt.outliers.channel_significance
            S.max_iter = 5;

            Dold=spm_eeg_load(spm_file);
            S2=[];
            S2.D=Dold;
            fname=fnamedat(Dold);
            [pth nm]=fileparts(fname);
            S2.outfile=[pth filesep 'S' nm '.mat'];
            D2=spm_eeg_copy(S2);
            D2 = osl_detect_artefacts(D2,S);
            D2.save();

            figs = [report.bad_channels(D2) report.bad_trials(D2)];
            opt_report=osl_report_set_figs(opt_report,arrayfun(@(x) sprintf('epoched_%s',get(x,'tag')),figs,'UniformOutput',false),figs,arrayfun(@(x) get(x,'name'),figs,'UniformOutput',false));
            opt_report=osl_report_print_figs(opt_report);

            % delete obsolete spm file
            spm_file_old=[opt.dirname filesep spm_files_epoched_basenames{subnum}];
            Dold=spm_eeg_load(spm_file_old);
            
            % in case we had loaded in an SPM file from the start, check it
            % is not the exact same name as old SPM file?
            
            if strcmp(opt.input_file_type,'spm_files')
            [p fname] = fileparts(opt.spm_files{subnum});
            if ((opt.cleanup_files == 1) || (opt.cleanup_files == 2)) && ~strcmp(Dold.fname, [fname '.mat'])
                Dold.delete;
            end;
            else
                % delete anyway?
                Dold.delete;
            end
                

            spm_files_epoched_basenames{subnum}=['S' spm_files_epoched_basenames{subnum}];

            opt_results.spm_files_epoched_basename=spm_files_epoched_basenames{subnum};

        end

        opt_results.spm_files_basename=spm_files_basenames{subnum};

        %%%%%%%%%%%%%%%%%%%
        %% generate source recon web report for this session
        opt_report=osl_report_write(opt_report);
        opt_report=osl_report_add_sub_report(opt_report, opt_report);

        %%%%%%%%%%%%%%%%%%%
        %% save opt results
        opt_results.fname=['session' num2str(subnum)];
        disp(['Saving opt results for: ' opt_results.fname]);
        opt_save_results(opt, opt_results);
        opt.results.fnames{subnum}=opt_results.fname;

    catch ME,

         % set output SPM file as empty for this subject:
        spm_files_basenames{subnum}=[];
        if(opt.epoch.do),
            spm_files_epoched_basenames{subnum}=[];
        end;

        disp(['Subject ' num2str(subnum) ' failed']);
        ME.getReport

    end
    
    % build logfile
    % the original command was 'cat opt_results.logfile '>' opt.results.logfile'
    % which implies that opt.results.logfile is overwritten because it uses > not >>
    % Thus replacing with a copyfile call here. Maybe a more sophisticated solution is
    % required if 
    copyfile(opt_results.logfile,opt.results.logfile);

end

diary(opt.results.logfile);

%%%%%%%%%%%%%%%%%%%%
%% gather results over all sessions
opt=opt_gather_results(opt);

%%%%%%%%%%%%%%%%%%%%
%% diagnostic plots over all sessions
[opt, opt_report]=opt_report_summary_plots(opt, opt_report);

%%%%%%%%%%%%%%%%%%%
%% generate web report
opt.results.report=osl_report_write(opt_report);

%%%%%%%%%%%%%%%%%%%
%% output res
opt.osl2_version=osl_version;

opt.fname=[opt.dirname filesep 'opt'];

save(opt.fname, 'opt');

disp(['To view OPT report, point your browser to <a href="' opt.results.report.html_fname '">' opt.results.report.html_fname '</a>']);

diary off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL EXTRAS BELOW:

%%%%%%%%%%%%%%%%%%%
%% Run this bit of code if you want to unreject all trials and channels
if(0)
    for subi= 1:length(spm_files), subnum=opt.sessions_to_do(subi);
        S2=[];
        spmfilename=spm_files{subnum};
        S2.D = spmfilename;
        D=spm_eeg_load(S2.D);

        rejected=sum(D.badtrials)
        tmp=1:length(D.conditions);
        trlsel = zeros(1, length(tmp));
        D = badtrials(D, tmp, trlsel);
        rejected=sum(D.badtrials)

        chans=D.indchantype('MEEG');
        D = D.badchannels(1:size(chans),0);
        D.save;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=set_bad(D,BadEpochs)

% Save bad epochs using method meeg/events
BadEvents = struct([]);
for ev = 1:numel(BadEpochs)
  if numel(BadEpochs{ev} == 2)
    BadEvents(ev).type     = 'artefact_OSL';
    BadEvents(ev).value    = 'all';
    BadEvents(ev).time     =  BadEpochs{ev}(1);
    BadEvents(ev).duration = diff(BadEpochs{ev});
    BadEvents(ev).offset = 0;
  end
end

% Load events
Events = D.events;

% Remove previous bad epoch events
if isfield(Events,'type')
  Events(strcmp({Events.type},'artefact_OSL')) = [];
end

% Concatenate new and old events
if size(Events,1) < size(Events,2)
  BadEvents = BadEvents(:);
end
if ~isempty(BadEvents)
  Events = [Events(:); BadEvents(:)];
end

% Save new events with previous
D = events(D,1,Events);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig_handles fig_names fig_titles] = plot_bad(D, S)

BadEpochs=S.BadEpochs;
fig_handles=[];
fig_names ={};
fig_titles={};

for ss=1:length(S.outlier_measure_fns),
    fig_handles(ss)=sfigure;
    set(fig_handles(ss),'Position',[1 1 1300 450]);
    for mm=1:length(S.modalities),
        chan_inds = find(strcmp(D.chantype, S.modalities{mm}));
        subplot(2,1,mm);
        tmpdat=squeeze(D(chan_inds,:,1));
        plot(D.time, feval(S.outlier_measure_fns{ss}, tmpdat));
        a=axis;
        ylims=a(3:4);

        for b = 1:numel(BadEpochs)
            line([BadEpochs{b}(1) BadEpochs{b}(1)],ylims,'linewidth',1.5,'linestyle','-','color',[0.1 0.8 0.1])
            if numel(BadEpochs{b}) == 2
              line([BadEpochs{b}(2) BadEpochs{b}(2)],ylims,'linewidth',1.5,'linestyle','-','color','r')
            end
        end

        plot4paper('time(s)',[S.outlier_measure_fns{ss} '(' S.modalities{mm} ')']);
    end;

    fig_names{ss}=['bad_segments_' S.outlier_measure_fns{ss} '_chans_with_bad_epochs_tc'];
    fig_titles{ss}=['BAD SEGMENTS: ' S.outlier_measure_fns{ss} ' over chans'];

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [maxfilter_failed D report]=opt_maxfilter_check_output(opt, fif_sss, spm_files_basename, report)

D=[];

maxfilter_failed=false;
if ~exist([fif_sss '.fif'],'file'),
    maxfilter_failed=true;
    return;
end;

%% Convert to SPM
disp(['%%%%%%%%%%%%%%%%%%%%%%%  CONVERT SSS-ed DATA, SESS = ' spm_files_basename '  %%%%%%%%%%%%%%%%%%%%%%%'])

S2=[];
spm_file=[opt.dirname filesep spm_files_basename];
S2.spm_file=spm_file;
[p fifname e] = fileparts(fif_sss);
S2.fif_file=[p filesep fifname '.fif'];

if(isfield(opt.convert,'trigger_channel_mask'))
    S2.trigger_channel_mask=opt.convert.trigger_channel_mask;
end;

[ D tmp fig_handles fig_names ] = osl_import(S2);

% test to make sure that all data is not zero-ed
chan_list=[];
for mm=1:length(opt.modalities),
    chan_list=[chan_list find(strcmp(chantype(D),opt.modalities{mm}))];
end;

if range(squash(D(chan_list,:,:)))==0
    maxfilter_failed=true;
end;

if maxfilter_failed
    return;
else

    if ~isempty(fig_handles),
        fig_titles={}; fig_titles{1}='CONVERT: histogram of trigger codes';
        report=osl_report_set_figs(report,fig_names,fig_handles,fig_titles);
        report=osl_report_print_figs(report);
    end;

    % delete fif file in opt dir that is no longer needed
    if(opt.cleanup_files == 2)
        delete(S2.fif_file);
    end;

end;

