function bashfilename = osl_prepare_oat_batch(oat,deploy,maxRam)
%
% function bashfilename = osl_prepare_oat_batch(oat,deploy,maxRam)
%
% Set deploy = 1 for mcc compiled code, deploy = 0 to run as matlab
% sessions (NB may use too many licences!)
%
% Set maxRam to an estimate of the maximum RAM each process will need, in
% Gb.  This will depend on the analysis.
%
% gw 2013

%% get osl going
osldir = '/home/xpsy-brain-cognition/gwallis/matlab/osl1.3.0';
addpath(osldir);
osl_startup(osldir);

%% load the oat
oat = osl_load_oat(oat);

%% work out which machine we're on and set the parameters for the job accordingly
res = [];
imachine = 0;
machines = {'hal','sal','arcus'};
while isempty(res)
    imachine = imachine + 1;
    res = runcmd(['hostname -f | grep '  machines{imachine}]);
end % while isempty(res)

if isempty(res)
    error('Can''t work out which machine this is running on!');
end

machine_node_ram        = [16,24,64]; % hal, sal, arcus
machine_node_ncores     = [8,8,16];   % hal, sal, arcus
nodeLimit               = [80,64,84];
ram                     = machine_node_ram(imachine);
ncores                  = machine_node_ncores(imachine);
coresPerNode            = min( ncores,floor(ram ./ maxRam) );
nSess                   = numel(oat.source_recon.sessions_to_do);
nNodes                  = uint8(ceil(nSess / coresPerNode));  % how many nodes are needed?
if nNodes > nodeLimit(imachine)
    error(['This job would need ' num2str(nNodes) ' nodes, but the machine only has ' num2str(nodeLimit(imachine)) ' nodes.']);
end % if nNodes > nodeLimit(imachine)

%% We can only run parallel matlab jobs WITHIN a node.  That means we need
% a bunch of shell scripts, each of which submits a limited number of jobs
% within the environs of a single node.  Then we launch these by calling a
% 'master script'

[oatdir,~] = fileparts(oat.fname);
[~,oatname] = fileparts(oatdir);

% open a 'launcher' file which just calls each of the batch scripts
launcherfilename = fullfile(oatdir,[oatname '_batcher.sh']);
launcherfile = fopen(launcherfilename,'w');
fprintf(launcherfile,'#!/bin/bash\n\n');

session_batch_lims1 = 1:nNodes;
session_batch_lims1 = session_batch_lims1 - 1;
session_batch_lims1 = 1 + session_batch_lims1*coresPerNode;
session_batch_lims2 = session_batch_lims1 + coresPerNode -1;
for iScript = 1:nNodes
    
    bashfilename = fullfile(oatdir,[oatname '_' num2str(iScript) 'of' num2str(nNodes) '.sh']);
    jobname = [oatname '_' num2str(iScript) 'of' num2str(nNodes)];
    bashfile = fopen(bashfilename,'w');
    
    % write the general stuff
    fprintf(bashfile,'#!/bin/bash\n\n');
    
    fprintf(bashfile,'#PBS -l nodes=1\n');
    fprintf(bashfile,'#PBS -l walltime=40:00:00\n');
    fprintf(bashfile,'#PBS -N %s\n\n',jobname);
    fprintf(bashfile,'#PBS -V \n\n');
    
    fprintf(bashfile,'# Go to the directory where the job was submitted\n');
    fprintf(bashfile,'cd $PBS_O_WORKDIR\n\n');
    
    fprintf(bashfile,'# Set up the FSL environment\n');
    fprintf(bashfile,'. $FSLDIR/etc/fslconf/fsl.sh\n\n');
    
    fprintf(bashfile,'# Make sure the executable is on the path\n');
    fprintf(bashfile,['export PATH=$HOME/bin:${PATH}\n\n']);
    
    for iS = session_batch_lims1(iScript) : session_batch_lims2(iScript)
        
        oattemp = oat;
        oattemp.to_do = [1 1 0 0];
        masterDir = oattemp.source_recon.dirname;
        [rootDir,master] = fileparts(masterDir);
        nam = [master '_sess' num2str(iS) '.oat'];
        subdir = fullfile(rootDir,nam);
        mkdir(subdir);
        oattemp.source_recon.dirname = subdir;
        oattemp.source_recon.sessions_to_do    = [iS];
        oattemp.first_level.sessions_to_do     = [iS];
        oattemp = osl_check_oat(oattemp);
        oattemp = osl_save_oat(oattemp);
        
        if deploy
            cmmd = ['runOSCoat "' oattemp.source_recon.dirname '" "' osldir '" &'];
        else
            cmmd = ['matlab -nodisplay -nosplash -r runOSCoat(''' oattemp.source_recon.dirname ''',''' osldir ''') &'];
        end
        fprintf(bashfile,'%s\n',cmmd);
    end % for iS = 1:nSess
    
    if deploy
        % keep the shell script open whilst the processes are running
        fprintf(bashfile,'\n');
        fprintf(bashfile,'# Keep the shell script open whilst the processes run\n');
        fprintf(bashfile,'while [ "$(pidof runOSCoat)" ]; do\n');
        fprintf(bashfile,'\tsleep 60\n');
        fprintf(bashfile,'done\n');
    else
        % keep the shell script open whilst the processes are running
        fprintf(bashfile,'\n');
        fprintf(bashfile,'# Keep the shell script open whilst the processes run\n');
        fprintf(bashfile,'while [ "$(pidof MATLAB)" ]; do\n');
        fprintf(bashfile,'\tsleep 60\n');
        fprintf(bashfile,'done\n');
    end
    
    fclose(bashfile);
    
    cmmd2 = ['qsub ' bashfile ' &'];
    fprintf(launcherfile,'%s\n',cmmd2);
    
end % for iScript = 1:nNodes

fclose(launcherfile);
disp(['The launcher file is at ' launcherfilename]);