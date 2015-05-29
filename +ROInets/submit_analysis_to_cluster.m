function [] = submit_analysis_to_cluster(Settings, DobjNames, logfileDir)
%SUBMIT_ANALYSIS_TO_CLUSTER  submit a set of sessions for network analysis
%
% SUBMIT_ANALYSIS_TO_CLUSTER(OIL) runs individual network analyses on
%   FMRIB cluster. 
%
% See also: osl_network_analysis


%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 373 $
%	$LastChangedDate: 2015-01-12 16:36:30 +0000 (Mon, 12 Jan 2015) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 01-May-2014 12:17:57


fprintf('%s: submitting analysis to cluster. \n', mfilename);

%% Parse input settings
fprintf('%s: checking inputs. \n', mfilename);

Settings  = ROInets.check_inputs(Settings);
nSessions = length(DobjNames);

% make results directory
ROInets.make_directory(Settings.outputDirectory);

% save settings
settingsFileName = fullfile(Settings.outputDirectory, 'ROInetworks_settings.mat');
save(settingsFileName, 'Settings');


%% Generate matlab script for each session
fprintf('Generating matlab scripts. \n');
for iSession = nSessions:-1:1,
    fprintf('Generating scripts to submit session %d out of %d. \n', nSessions-iSession+1, nSessions);

    m_scriptName{iSession} = fullfile(Settings.outputDirectory,                                ...
                                      sprintf('matlab_submit_individual_network_session_%d.m', ...
                                              iSession));
    fid = fopen(m_scriptName{iSession}, 'w');
    if -1 == fid, 
        error([mfilename ':fileWritingFailed'], ...
              'Could not open %s for writing. \n', ...
              m_scriptName{iSession});
    else
        % neat file cleanup - will close the file, even if the script
        % errors
        C = onCleanup(@() fclose(fid));
    end%if
    
    % load D object to be analysed
    Dfilt = spm_eeg_load(DobjNames{iSession});
    Dfilt = Dfilt.montage('switch', 2);
    Dfilt.save;
    
    % set a save file name
    sessionName = strrep(Dfilt.fname, '.mat', '_ROInets');
    resultsName = fullfile(Settings.outputDirectory, sessionName);
    
    % make file header
    fprintf(fid, '%% matlab_submit_individual_network_session_%d.m\n', iSession);
    fprintf(fid, '%% Submits network analysis for session %d\n', iSession);
    fprintf(fid, '%% Made on %s\n', datestr(now));
    fprintf(fid, '\n\n');
    
    % check for existing results - we're not going to over-write.
    fprintf(fid, '%% Check to see if results already computed. \n');
    fprintf(fid, 'if exist(''%s.mat'', ''file''), exit; end\n\n', fullfile(Settings.outputDirectory, sessionName));
    
    % set up paths and osl
    fprintf(fid, '%% Open OSL\n');
    fprintf(fid, 'set(0, ''DefaultFigureVisible'', ''off'');\n');
    fprintf(fid, 'OSLDIR = ''/home/fs0/gilesc/OSL-Repo/osl2-outer/'';\n');
    fprintf(fid, 'addpath(OSLDIR);\n');
    fprintf(fid, 'osl_startup(OSLDIR);\n\n');
    
    % Giles' paths
    fprintf(fid, '%% Set up paths\n');
    fprintf(fid, 'addpath(genpath(''/home/fs0/gilesc/IBME-SVN-Repo/src/''));\n');
    
    % Load information
    fprintf(fid, '%% Load data for analysis\n');
    fprintf(fid, 'fprintf(''Loading data for analysis \\n'');\n');
    fprintf(fid, 'tmp      = load(''%s'');\nSettings = tmp.Settings;\n', settingsFileName);
    fprintf(fid, 'clear tmp\n');
    fprintf(fid, 'Settings.sessionName = ''%s'';\n', sessionName);
    
    % Run analysis
    fprintf(fid, '\n%% Run the analysis\n');
    fprintf(fid, 'fprintf(''Running run_individual_correlation_analysis \\n'');\n');
    fprintf(fid, 'mats = ROInets.run_individual_network_analysis(''%s'', Settings, ''%s'');\n', Dfilt.fnamedat, resultsName);
    fprintf(fid, 'fprintf(''All done. \\n'');\n');
    % close matlab session
    fprintf(fid, '\nexit\n\n');
    
    % close file
    clear C
end%loop over sessions

%% Generate submit script for each session
fprintf('Generating individual submit scripts. \n');
for iSession = nSessions:-1:1,
    fprintf('Generating scripts to submit session %d out of %d. \n', nSessions-iSession+1, nSessions);

    bash_scriptName{iSession} = fullfile(Settings.outputDirectory,                                ...
                                    sprintf('shell_submit_individual_network_session_%d', ...
                                            iSession));
    fid = fopen(bash_scriptName{iSession}, 'w');
    if -1 == fid, 
        error([mfilename ':fileWritingFailed'], ...
              'Could not open %s for writing. \n', ...
              bash_scriptName{iSession});
    else
        % neat file cleanup - will close the file, even if the script
        % errors
        C = onCleanup(@() fclose(fid));
    end%if
    
    % file header
    fprintf(fid, '#! /bin/sh\n');
    fprintf(fid, '#%s\n', bash_scriptName{iSession});
    fprintf(fid, '# Script to run %s on the cluster\n', m_scriptName{iSession});
    fprintf(fid, '# Made on %s\n', datestr(now));
    fprintf(fid, '\n\n');
    
    % deal with virtual display
    fprintf(fid, '/usr/local/bin/virtual-x -f > /tmp/$$.virtual-display\n');
    fprintf(fid, 'disp=`grep ''display started'' /tmp/$$.virtual-display | awk -F: ''{print $NF}''`\n');
    fprintf(fid, 'rm /tmp/$$.virtual-display\n');
    fprintf(fid, 'export DISPLAY=":$disp"\n');
    fprintf(fid, '\n');
    logFileName = fullfile(Settings.outputDirectory, sprintf('shell_submit_individual_network_session_%d_logfile', iSession));
    fprintf(fid, 'matlab -singleCompThread -nodesktop -logfile %s -nosplash \\< %s\n', logFileName, m_scriptName{iSession});
    fprintf(fid, '\n');
    fprintf(fid, 'virtual-x -k $disp\n');
    
    % close file
    clear C
end%loop over sessions
 

%% Submit jobs using one master script
fprintf('Generating master script. \n');
submitFileName = fullfile(Settings.outputDirectory, ...
                          'submit_all_individual_network_analyses');
                      
fid = fopen(submitFileName, 'w');
if -1 == fid,
    error([mfilename ':fileWritingFailed'], ...
          'Could not open %s for writing. \n', ...
          submitFileName);
else
    % neat file cleanup - will close the file, even if the script
    % errors
    C = onCleanup(@() fclose(fid));
end%if

% header
fprintf(fid, '#! /bin/sh\n');
fprintf(fid, '# %s\n', submitFileName);
fprintf(fid, '# Script to submit all individual network analyses to cluster\n');
fprintf(fid, '# Made on %s\n', datestr(now));
fprintf(fid, '\n\n');

% submit commands
for iSession = 1:nSessions,
    jobName = sprintf('ROInetworks-sub%d', iSession);
    logName = fullfile(logfileDir, sprintf('fsl_networks_sub_log_sess_%d', iSession));
    % make file executable
    fprintf(fid, 'chmod a+x %s\n', bash_scriptName{iSession});
    % submit job to queue
    fprintf(fid, 'fsl_sub -q bigmem.q -N %s -l %s %s\n', jobName, logName, bash_scriptName{iSession});
end%loop over sessions

% close file
clear C

% run the file
fprintf('Submitting ..... \n');
GC_checked_system(sprintf('chmod a+x %s', submitFileName));
% GC_checked_system(sprintf('source %s',    submitFileName))

fprintf('   Done. \n');

end%submit_analysis_to_cluster
