function oat = osl_run_oat(oat)

% oat = osl_run_oat(oat)
%
% Runs an OAT (OHBA's easy Analysis Tool) using the settings contained in
% oat, which needs to be setup by calling oat=osl_setup_oat(S);
%
% MWW 2011

OSLDIR = getenv('OSLDIR');

% check settings
oat = osl_check_oat(oat);

% make dir for the results:
if(isempty(findstr(oat.source_recon.dirname, '.oat')))
    oat.source_recon.dirname=[oat.source_recon.dirname, '.oat'];
end

if ~exist(oat.source_recon.dirname,'dir')
    mkdir(oat.source_recon.dirname);
end

% set results container
oat.results=[];
oat.results.date=date;

% set logfile up
oat.results.plotsdir=fullfile(oat.source_recon.dirname,'plots');
oat.results.logfile=fullfile(oat.results.plotsdir,['log-' date '.txt']);
% delete any existing diary file with the same name
delete(oat.results.logfile);

% set first level diagnostic report up
oat.results.report=osl_report_setup(oat.results.plotsdir,['OAT report'],oat.results.logfile);

diary(oat.results.logfile);

% turns figure plotting on or off in SPM calls
spm_get_defaults('cmdline',~oat.do_plots);

try
    if(oat.to_do(1))

        disp('*************************************************************');
        disp('Running source_recon');
        disp('*************************************************************');

        
        if(strcmp(oat.source_recon.method,'none'))
            [results_fnames results]=oat_run_source_recon_sensorspace(oat);
        else
            [results_fnames results]=oat_run_source_recon(oat);
%         else
%             error('recon_method is invalid');
        end

        % for those that have just been run, overwrite results_fnames:
        oat.source_recon.results_fnames(oat.source_recon.sessions_to_do)=results_fnames(oat.source_recon.sessions_to_do);

        oat.results.report=osl_report_add_sub_report(oat.results.report, results.report);

        % save OAT to disk
        oat = osl_save_oat(oat);

    end
catch ME
    disp('Source recon OAT stage has failed');
    ME.getReport
    return;
end

%try,
    if(oat.to_do(2))
    
        disp('*************************************************************');
        disp('Running first_level');
        disp('*************************************************************');

        try
            tmp=oat.source_recon.results_fnames;
        catch
            try
                warning('oat does not containing any results for previous level. Trying to load in results using oat=osl_load_oat(oat).');
                oatin=osl_load_oat(oat);
                oat.source_recon.results_fnames=oatin.source_recon.results_fnames;
                disp(['Loaded in ' oatin.fname]);
            catch
                error('Failed: Source recon stage has not been run.');
            end
        end

        [results_fnames results] = oat_run_first_level(oat);

        % for those that have just been run, overwrite results_fnames:
        oat.first_level.results_fnames(oat.first_level.sessions_to_do)=results_fnames(oat.first_level.sessions_to_do);
        
        oat.results.report=osl_report_add_sub_report(oat.results.report, results.report);

        % save OAT to disk
        oat = osl_save_oat(oat);

    end
%catch ME,
%    disp('First level OAT stage has failed');
%    ME.getReport
%    return;
%end;

try
    if(oat.to_do(3))

        disp('*************************************************************');
        disp('Running subject_level');
        disp('*************************************************************');

        try
            tmp = oat.first_level.results_fnames;
        catch
            try
                warning('oat does not containing any results for previous levels. Trying to load in results using oat=osl_load_oat(oat).');
                oatin=osl_load_oat(oat);
                oat.first_level.results_fnames=oatin.first_level.results_fnames;
                disp(['Loaded in ' oatin.fname]);
            catch
                error('Failed: First level stage has not been run.');
            end

        end

        results_fnames = oat_run_subject_level(oat);

        oat.subject_level.results_fnames(oat.subject_level.subjects_to_do) = results_fnames(oat.subject_level.subjects_to_do);

        % note there is no report for subject level
        
        % save OAT to disk
        oat = osl_save_oat(oat);
    end
catch ME
    disp('Subject level OAT stage has failed');
    ME.getReport
    return;
end

try
    if(oat.to_do(4))

        disp('*************************************************************');
        disp('Running group_level');
        disp('*************************************************************');

        try
            tmp = oat.subject_level.results_fnames;
        catch
            try
                warning('oat does not containing any results for previous level. Trying to load in results using oat=osl_load_oat(oat).');
                oatin=osl_load_oat(oat);
                oat.subject_level.results_fnames=oatin.subject_level.results_fnames;
                disp(['Loaded in ' oatin.fname]);
            catch
                error('Failed: Subject level stage has not been run.');
            end
        end

        [oat.group_level.results_fnames results]= oat_run_group_level(oat);

        oat.results.report=osl_report_add_sub_report(oat.results.report, results.report);
     
        disp(['To view group level report, point your browser to <a href="' results.report.html_fname '">' results.report.html_fname '</a>']);

        % save OAT to disk
        oat = osl_save_oat(oat);
    end
catch ME
    disp('Group level OAT stage has failed');
    ME.getReport
    return;
end

%%%%%%%%%%%%%%%%%%%
%% generate web report
oat.results.report=osl_report_write(oat.results.report);
disp(['To view OAT report, point your browser to <a href="' oat.results.report.html_fname '">' oat.results.report.html_fname '</a>']);

diary off;
