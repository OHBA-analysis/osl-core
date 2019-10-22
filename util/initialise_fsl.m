function initialise_fsl()
% Adds FSL directories to system path, sets FSL environment variables
% and adds Matlab utilities to the path
%
% Reads FSL directories from the OSL configuration file to provide support for multiple
% versions
% 
% Specifically, it uses the following variables in osl.conf:
%
% FSLDIR
% FSLBIN
% FSLLIB
    
    s = osl_conf.read();

    % Initialize FSL
    setenv('FSLDIR',s.FSLDIR)
    setenv('FSLOUTPUTTYPE','NIFTI_GZ')
    if ~osl_util.contains(getenv('PATH'),s.FSLBIN)
        setenv('PATH',sprintf('%s%s%s',s.FSLBIN,pathsep,getenv('PATH')));
    end
    if ~osl_util.contains(getenv('LD_LIBRARY_PATH'),s.FSLLIB)
        setenv('LD_LIBRARY_PATH',sprintf('%s%s%s',s.FSLLIB,pathsep,getenv('LD_LIBRARY_PATH')));
    end

    % Add FSLDIR to path, if it exists - otherwise, print warning
    if osl_util.isdir(getenv('FSLDIR'))
        addpath(fullfile(s.FSLDIR,'etc','matlab'));
    else
        fprintf(2,'FSLDIR does not exist. Check that it is set correctly in osl.conf\n');
    end
    
    if isempty(which('read_avw'))
        fprintf(2,'FSL matlab utilities not found\n');
    end
    
    [status,res] = system('fslval');
    if status~=1
        fprintf(2,'FSL is not installed properly. Have you installed the correct version and set it in osl.conf?\n');
    end

end
