function osl_run_parallel_oat(oat,osldir)
%
% function runOSCoat(oat,osldir)
%
% Loads and runs a prepared oat analysis.  Compile with osl_make_standalone
% to use on the OSC.

disp('Check input:');
class(oat)
class(osldir)

try
    
    %% set up the paths and initialize osl
    if isdeployed
        % the path to osl_dir should already be included as part of the
        % compilation
        osl_startup(osldir);
    else
        addpath(osldir);
        osl_startup(osldir);
    end
    oat = osl_load_oat(oat);
    sess_num = oat.source_recon.sessions_to_do;
    logfilename = [oat.source_recon.dirname '/log_sess_' num2str(sess_num) '.txt'];
    diary(logfilename)
    disp(['About to run oat for session number ' num2str(sess_num)]);
    osl_run_oat(oat);
    diary off
    
catch caught_error
    
    if isstruct(oat)
        fn = (fullfile(oat.source_recon.dirname,'caught_error'));
        save(fn,'caught_error');
    else
        fn = (fullfile(oat,'caught_error'));
        save(fn,'caught_error');
    end
    
end
