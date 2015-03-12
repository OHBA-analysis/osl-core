function [ command ] = cleanEnvironment(  )
% Returns command to clean up the environment variables MATLAB shits all
% over before running an external command
%
% Sam Harrison

command = ['unset  XKEYSYMDB TERM LC_NUMERIC XAPPLRESDIR AUTOMOUNT_MAP' ...
        ' GFORTRAN_STDIN_UNIT OSG_LD_LIBRARY_PATH TOOLBOX LC_MESSAGES' ...
        ' XFILESEARCHPATH BASEMATLABPATH DYLD_FRAMEWORK_PATH' ...
        ' DYLD_LIBRARY_PATH GFORTRAN_STDERR_UNIT ARCH; '];

end

