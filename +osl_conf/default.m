function s = default()
%
% Default values for required config fields.
%   
% JH

    [fsldir,fslbin,fsllib] = default_fsldir();
    
    s.FSLDIR = fsldir;
    s.FSLBIN = fslbin;
    s.FSLLIB = fsllib;
    
    
    s.SPMDIR = '';
    s.WORKBENCH = '';

end

function [fsldir,bindir,libdir] = default_fsldir()
% Guess where the FSL directories maybe

    fsldir = getenv('FSLDIR');
    
    if ~isempty(fsldir) % manually specified
        bindir = fullfile(fsldir,'bin');
        libdir = fullfile(fsldir,'lib');
        
    elseif osl_util.isdir('/usr/local/fsl') % default location
        fsldir = '/usr/local/fsl';
        bindir = '/usr/local/fsl/bin';
        libdir = '/usr/local/fsl/lib';
        
    elseif osl_util.isdir('/usr/share/fsl/5.0') % Debian FSL package
        fsldir = '/usr/share/fsl/5.0';
        bindir = '/usr/share/fsl/5.0/bin';
        libdir = '/usr/lib/fsl/5.0';
        
    else
        warning( 'Could not find FSLDIR.' );
    end
    
end