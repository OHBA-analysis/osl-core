function osl_startup( osldir )

global OSLDIR;

% osl_startup( osldir )
%
% SETS UP THE BASIC PATH SETTINGS
%
% MWW 2012
%
% does no path-changing if running in deployed mode (gw '13).

if ~isdeployed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check path for other versions of:
    checklist={'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt'};
    
    % and remove them:
    oldpath=path;
    indcolon=findstr(oldpath,':');
    st=1;
    restoredefaultpath;
    
    restoredpath=path;
    addpath(genpath(osldir));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check fsl has been setup
    if isempty(getenv('FSLDIR')),
       error('The environmental variable FSLDIR is not set. Please exit Matlab, ensure FSLDIR is set and then restart Matlab. See the Prerequisites section at https://sites.google.com/site/ohbaosl/download'); 
    end;

    % try a dummy call to an fsl tool to make sure FSL is properly
    % installed:
    [status,res]=dos(['fslval']);
    if(status~=1)
        error('FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab. See the Prerequisites section at https://sites.google.com/site/ohbaosl/download');
    end;
    
    addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    jj=1;
    found=zeros(length(checklist),1);
    
    for ii=1:length(indcolon),
        pathlistitem=oldpath(st:indcolon(ii)-1);
        
        ok=1;
        for kk=1:length(checklist),
            if(any(findstr(pathlistitem,checklist{kk}))),
                ok=0;
                found(kk)=1;
            end;
        end;
        
        if(ok)
            if(~any(findstr(restoredpath,pathlistitem)))
                addpath(pathlistitem);
            end;
        end;
        
        st=indcolon(ii)+1;
    end;
    
end % if ~ isdeployed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copy changes to SPM code from osl
filelist={};targetdir={};

filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_output_montage_osl.m';
targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_write_spmeeg_osl.m';
targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

filelist{end+1}='osl2/spm-changes/ft_read_event_4osl.m';
targetdir{end+1}='spm12/external/fieldtrip/fileio';

filelist{end+1}='osl2/spm-changes/read_trigger_4osl.m';
targetdir{end+1}='spm12/external/fieldtrip/fileio/private';

filelist{end+1}='osl2/spm-changes/badsamples.m';
targetdir{end+1}='spm12/@meeg';

% filelist{end+1}='osl2/spm-changes/lmoutr.mexmaci64';
% targetdir{end+1}='spm12/external/fieldtrip/utilities/private';
% 
% filelist{end+1}='osl2/spm-changes/mxSerialize.mexmaci64';
% targetdir{end+1}='spm12/external/fieldtrip/utilities/private';
% 
% filelist{end+1}='osl2/spm-changes/ptriproj.mexmaci64';
% targetdir{end+1}='spm12/external/fieldtrip/utilities/private';

filelist{end+1} = 'osl2/spm-changes/ft_headmodel_localspheres.m';
targetdir{end+1}='spm12/external/fieldtrip/forward/';

filelist{end+1}='osl2/spm-changes/path.m';
targetdir{end+1}='spm12/@meeg';

for kk=1:length(filelist),
    runcmd(['cp -f ' osldir '/' filelist{kk} ' ' osldir '/' targetdir{kk}]);
    rmpath([osldir '/' targetdir{kk}]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OSLDIR=osldir;

spm_get_defaults('cmdline',true);
spm eeg;
close all;

if ~isdeployed
    for kk=1:length(checklist),
        if found(kk),
            warning(['Found and removed paths that contain the string ' checklist{kk}]);
        end;
    end;
end % if ~ isdeployed
