function results = osl_run_maxfilter(settings_in)
     
% results = osl_run_maxfilter(settings)

error('This function is being removed, and will not work if run');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% required inputs
    
try, settings.raw_fif_files=settings_in.raw_fif_files; settings_in = rmfield(settings_in,'raw_fif_files'); catch, settings.raw_fif_files=[]; end; % Specify a list of the raw fif files for subjects

% check list of SPM MEEG filenames input
if(~isempty(settings.raw_fif_files))
    sess=settings.raw_fif_files;    
end

num_sessions=length(sess);

% check that full directory names have been specified
for iSess = 1:num_sessions,
    sessPath = fileparts(sess{iSess});
    if isempty(sessPath) || strcmpi(sessPath(1), '.'),
        error([mfilename ':FullPathNotSpecified'], ...
              'Please specify full paths for the fif, input or spm files. \n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optional settings:

try, settings.sessions_to_do=settings_in.sessions_to_do; settings_in = rmfield(settings_in,'sessions_to_do'); catch, settings.sessions_to_do=1:num_sessions; end;

try, settings.dirname=settings_in.dirname; settings_in = rmfield(settings_in,'dirname'); catch, settings.dirname=[sess{1} '.maxfilter']; end; % directory name 
if(isempty(findstr(settings.dirname, '.maxfilter')))
    settings.dirname=[settings.dirname, '.maxfilter'];
end

try, settings.remote_port=settings_in.maxfilter.remote_port; settings_in.maxfilter = rmfield(settings_in.maxfilter,'remote_port');
catch
    settings.remote_port = 0;
end

try, settings.do_sss=settings_in.maxfilter.do_sss; settings_in.maxfilter = rmfield(settings_in.maxfilter,'do_sss'); catch, settings.do_sss=1; end; % flag to indicate whether actual SSS maxfiltering should be done or not
try, settings.do_remove_badchans_pre_sss=settings_in.maxfilter.do_remove_badchans_pre_sss; settings_in.maxfilter = rmfield(settings_in.maxfilter,'do_remove_badchans_pre_sss'); catch, settings.do_remove_badchans_pre_sss=1; end; % flag to indicate whether bad chans should be removed before running SSS
try, settings.max_badchans_pre_sss=settings_in.maxfilter.max_badchans_pre_sss; settings_in.maxfilter = rmfield(settings_in.maxfilter,'max_badchans_pre_sss'); catch, settings.max_badchans_pre_sss=10; end; % maximum number of bad chans to be removed before running SSS
try, settings.movement_compensation=settings_in.maxfilter.movement_compensation; settings_in.maxfilter = rmfield(settings_in.maxfilter,'movement_compensation'); catch, settings.movement_compensation=1; end; % flag to indicate whether move comp should be done
try, settings.trans_ref_file=settings_in.maxfilter.trans_ref_file; settings_in.maxfilter = rmfield(settings_in.maxfilter,'trans_ref_file'); catch, settings.trans_ref_file=[]; end; % trans reference file to pass to maxfilter call using the -trans flag
try, settings.temporal_extension=settings_in.maxfilter.temporal_extension; settings_in.maxfilter = rmfield(settings_in.maxfilter,'temporal_extension'); catch, settings.temporal_extension=0; end; % flag to indicate whether Maxfilter temporal extension should be done
try, settings.maxfilt_dir=settings_in.maxfilter.maxfilt_dir; settings_in.maxfilter = rmfield(settings_in.maxfilter,'maxfilt_dir'); catch, settings.maxfilt_dir='/neuro/bin/util'; end; % where to find MaxFilter exe. Defaults to S.maxfilt_dir = '/neuro/bin/util'.
try, settings.bad_epochs=settings_in.maxfilter.bad_epochs; settings_in.maxfilter = rmfield(settings_in.maxfilter,'bad_epochs'); catch, settings.bad_epochs=cell(num_sessions,1); end; % Bad epochs to ignore (by maxfilter (passed using the -skip Maxfilter settingsion), one cell for each session, where the cell contains a (N_epochs x 2) matrix of epochs, where each row indicates the start and end time of each bad epoch (in secs)
try, settings.cal_file = settings_in.maxfilter.cal_file; settings_in.maxfilter = rmfield(settings_in.maxfilter,'cal_file'); catch, settings.cal_file = 0;end
try, settings.ctc_file = settings_in.maxfilter.ctc_file; settings_in.maxfilter = rmfield(settings_in.maxfilter,'ctc_file'); catch, settings.ctc_file = 0;end
try, settings.downsample.do=settings_in.downsample.do; settings_in.downsample = rmfield(settings_in.downsample,'do'); catch, settings.downsample.do=1; end; % flag to do or not do downsample
try, settings.downsample.freq=settings_in.downsample.freq; settings_in.downsample = rmfield(settings_in.downsample,'freq'); catch, settings.downsample.freq=250; end; % downsample freq in Hz


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['%%%%%%%%%%%%%%%%%%%%%%%  STARTING LOOP OVER SESSIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

results=[];
settings_report=[];
    
for subi=1:length(settings.sessions_to_do),
    subnum=settings.sessions_to_do(subi);

    % single session results container:

    %if(~settings.do_sss_maxfilter ||  (settings.do_remove_badchans_pre_sss && settings.do_sss_maxfilter))

    %%%%%%%%%%%%%%%%%%
    %% Create 1st round of non-SSSed data if doing double Maxfilter
    % OR create final non-SSSed Maxfiltered data if not doing SSS
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  MAXFILT WITH NO SSS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
    Smf=[];
    Smf.maxfilt_dir=settings.maxfilt_dir;
    [p fifname e] = fileparts(settings.raw_fif_files{subnum});
    Smf.fif=[p filesep fifname];

    if settings.ctc_file ~= 0
        Smf.ctc_file = ctc_file;
    end
    if settings.cal_file ~= 0
        Smf.cal_file = cal_file;
    end

    % output fif file from maxfilter call will be placed in the output
    % directory with a name based on settings.files_basenames{subnum}
    fifname_out=['nosss_fif_sesssion' num2str(subnum)];
    Smf.fif_out=[settings.dirname filesep fifname_out];
    Smf.logfile=1;

    if ~osl_util.isfile([Smf.fif '.fif']),
        error(['settings.raw_fif_files{' num2str(subnum) '} file ' [Smf.fif '.fif'] ' does not exist']);
    end

    Smf.movement_compensation=0;

    if(~Smf.movement_compensation)
        disp('Assuming data was sampled at 1000 Hz');
        sample_freq=1000;
        if(settings.downsample.freq<=sample_freq/4 && settings.downsample.do)
            Smf.downsample_factor=4;
        elseif(settings.downsample.freq<=sample_freq/2 && settings.downsample.do)
            Smf.downsample_factor=2;
        else
            Smf.downsample_factor=1;
        end
    end

    Smf.bad_epochs=settings.bad_epochs{subnum};

    Smf.nosss=1;
    if settings.remote_port ~= 0
        fif_sss=osl_call_maxfilter_remote(Smf,settings.remote_port);
    else
        fif_sss=osl_call_maxfilter(Smf);
    end

    type([fif_sss '_log.txt']);

    if ~osl_util.isfile([fif_sss '.fif'])
        warning([fif_sss '.fif does not exist. Can not continue with this subject.']);
    end

    %%%%%%%%%%%%%%%%%%%
    %% Convert to SPM

    S2=[];
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  CONVERT, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])

    spm_file=fullfile(settings.dirname,['nosss' settings.files_basenames{subnum}]);
    S2.spm_file=spm_file;
    [p fifname e] = fileparts(Smf.fif_out);
    S2.fif_file=[p filesep fifname '.fif'];
    if(isfield(settings.convert,'trigger_channel_mask'))
        S2.trigger_channel_mask=settings.convert.trigger_channel_mask;
    end
    [D] = osl_import(S2);

    % delete fif file in settings dir that is no longer needed
    if(settings.cleanup_files == 1) || (settings.cleanup_files == 2)
        delete(S2.fif_file);
    end

    close all

    %%%%%%%%%%%%%%%%%%%
    %%
    if(settings.do_sss)

        if(settings.do_remove_badchans_pre_sss),

            %%%%%%%%%%%%%%%%%%%
            %% Do DOUBLE Maxfilter procedure:
            %% First detects and remove bad channels on non-SSS data,
            %% then runs SSS
            disp(['%%%%%%%%%%%%%%%%%%%%%%%  RM BAD CHANS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
            spm_file=fullfile(settings.dirname ,['nosss' settings.files_basenames{subnum}]);
            D=spm_eeg_load(spm_file);

            for ii=1:length(settings.modalities),
                printprefix_mod=[settings.modalities{ii} '_preSSS'];

                %%%%%%%%%%%%%%%%%%%%%
                % use outlier detection approach
                if(1),

                    S = [];
                    S.D = D;
                    S.modalities{1}=settings.modalities{ii};
                    S.do_plot=0;
                    S.max_iter=1;
                    %S.outlier_measure_fns={'min','std'};
                    S.outlier_measure_fns={'std'};
                    S.outlier_wthresh_chan=[0.01];
                    S.just_chans=1;
                    S.plot_basename='double_maxfilter';
                    S.plot_basetitle='MAXFILTER: ';
                    S.max_bad_channels=settings.max_badchans_pre_sss;
                    [D2]=osl_detect_badevent(S);

                    list_bch=badchannels(D2);

                    if(~isempty(list_bch))
                        bad_list=badchannels(D);
                        vect_bch=zeros(1,size(D,1));
                        vect_bch(bad_list)=1;
                        vect_bch(list_bch)=1;

                        D = D.badchannels(cat(2,bad_list,list_bch),1);
                        D.save;
                    end

                    % delete obsolete spm file
                    spm_file_old=[settings.dirname filesep 'Snosss' settings.files_basenames{subnum}];
                    Dold=spm_eeg_load(spm_file_old);

                    if(settings.cleanup_files == 1) || (settings.cleanup_files == 2)
                        Dold.delete;
                    end

                    %%%%%%%%%%%%%
                    % diagnostic plots
                    S=[];
                    S.D=D;
                    S.print_plots=1;
                    S.plotsdir=settings.results.plotsdir;
                    S.modality=settings.modalities{ii};
                    S.printprefix=[printprefix_mod '_outlier'];
                    S.plot_basename='double_maxfilter';
                    S.plot_basetitle='MAXFILTER: ';
                    [res fig_names fig_handles fig_titles]=osl_check_bad_chans(S); % TODO - remove this call 

                    %settings_report=osl_report_set_figs(settings_report,fig_names,fig_handles,fig_titles);
                    %settings_report=osl_report_print_figs(settings_report);

                    %%%%%%%%%%%%%%%
                end
                %%%%%%%%%%%%%%%%%%%%%

                D.save;

            end

            results.badchans_pre_sss=badchannels(D);

        end

        %%%%%%%%%%%%%%%%%%%
        %% Re-Maxfilter data using bad channel information

        disp(['%%%%%%%%%%%%%%%%%%%%%%%  MAXFILT WITH SSS, SESS = ' num2str(subnum) '  %%%%%%%%%%%%%%%%%%%%%%%'])
        Smf=[];
        [p fifname e] = fileparts(settings.raw_fif_files{subnum});
        Smf.fif=[p filesep fifname];

        if settings.ctc_file ~= 0
            Smf.ctc_file = ctc_file;
        end
        if settings.cal_file ~= 0
            Smf.cal_file = cal_file;
        end

        % output fif file from maxfilter call will be placed in the settings
        % directory with a name based on settings.files_basenames{subnum}
        fifname_out=['sss_fif_' settings.files_basenames{subnum}];
        Smf.fif_out=[settings.dirname filesep fifname_out];

        Smf.logfile=1;
        Smf.nosss=0;

        spm_file=[settings.dirname filesep 'nosss' settings.files_basenames{subnum}]
        Smf.spmfile=spm_file;

        if ~isempty(settings.trans_ref_file) && ~strcmp(settings.trans_ref_file,''),
            [pth fifname ext] = fileparts(settings.trans_ref_file);
            Smf.trans_ref_file=[pth filesep fifname];
        end

        Smf.movement_compensation=settings.movement_compensation; %can only run with sss AND downsampling has to be switched off

        if(settings.downsample.do && Smf.movement_compensation)
            disp('Can not do any downsampling during Maxfilter call when Smf.movement_compensation is on.');
        end

        Smf.st=[];
        Smf.st.do=settings.temporal_extension;

        if(~Smf.movement_compensation)
            disp('Assuming data was sampled at 1000 Hz');
            sample_freq=1000;
            if(settings.downsample.freq<=sample_freq/4 && settings.downsample.do)
                Smf.downsample_factor=4;
            elseif(settings.downsample.freq<=sample_freq/2 && settings.downsample.do)
                Smf.downsample_factor=2;
            else
                Smf.downsample_factor=1;
            end
        end

        Smf.bad_epochs=settings.bad_epochs{subnum};

        if settings.remote_port ~= 0
            fif_sss=osl_call_maxfilter_remote(Smf,settings.remote_port);
        else
            fif_sss=osl_call_maxfilter(Smf);
        end
        type([fif_sss '_log.txt']);

        results.sss_autobad_off(subnum)=0;

        % check if Maxfilter has worked
        spm_files_basenames{subnum}=['sss' settings.files_basenames{subnum}];
        [maxfilter_failed D settings_report]=maxfilter_check_output(settings, fif_sss, spm_files_basenames{subnum}, settings_report);

        if maxfilter_failed,
            warning('Maxfilter has FAILED. Will retry with no bad channels set and autobad off!!! Results may be suspect.');

            results.sss_autobad_off(subnum)=1;
            % set all chans good
            D=spm_eeg_load(spm_file);
            chans=D.indchantype('MEEG');
            D = D.badchannels(1:size(chans),0);
            D.save;

            Smf.autobad_off=1;
            if settings.remote_port ~= 0
                fif_sss=osl_call_maxfilter_remote(Smf,settings.remote_port);
            else
                fif_sss=osl_call_maxfilter(Smf);
            end
            type([fif_sss '_log.txt']);

            % check if Maxfilter has worked
            [maxfilter_failed D settings_report]=maxfilter_check_output(settings, fif_sss, spm_files_basenames{subnum}, settings_report);

            if maxfilter_failed,
                error('Maxfilter has still FAILED. Can not continue with this subject.');
            end
        end

        %%%
        % delete obsolete nosss spm file with pre sss bad channels
        spm_file_old=[settings.dirname filesep 'nosss' settings.files_basenames{subnum}];
        Dold=spm_eeg_load(spm_file_old);
        if(settings.cleanup_files == 1) || (settings.cleanup_files == 2)
            Dold.delete;
        end
        %%%

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [maxfilter_failed D report]=opt_maxfilter_check_output(opt, fif_sss, spm_files_basename, report)

D=[];

maxfilter_failed=false;
if ~osl_util.isfile([fif_sss '.fif']),
    maxfilter_failed=true;
    return;
end

%% Convert to SPM
disp(['%%%%%%%%%%%%%%%%%%%%%%%  CONVERT SSS-ed DATA, SESS = ' spm_files_basename '  %%%%%%%%%%%%%%%%%%%%%%%'])

S2=[];
spm_file=[opt.dirname filesep spm_files_basename];
S2.spm_file=spm_file;
[p fifname e] = fileparts(fif_sss);
S2.fif_file=[p filesep fifname '.fif'];

if(isfield(opt.convert,'trigger_channel_mask'))
    S2.trigger_channel_mask=opt.convert.trigger_channel_mask;
end

[ D tmp fig_handles fig_names ] = osl_import(S2);

% test to make sure that all data is not zero-ed
chan_list=[];
for mm=1:length(opt.modalities),
    chan_list=[chan_list find(strcmp(chantype(D),opt.modalities{mm}))];
end

if range(squash(D(chan_list,:,:)))==0
    maxfilter_failed=true;
end

if maxfilter_failed
    return;
else

    if ~isempty(fig_handles),
        fig_titles={}; fig_titles{1}='CONVERT: histogram of trigger codes';
        %report=osl_report_set_figs(report,fig_names,fig_handles,fig_titles);
        %report=osl_report_print_figs(report);
    end

    % delete fif file in opt dir that is no longer needed
    if(opt.cleanup_files == 2)
        delete(S2.fif_file);
    end

end
