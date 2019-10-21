function initialise_workbench()
% Add workbench to system path
    
    s = osl_conf.read();
    
    if isempty(s.WORKBENCH), return; end
    
    wbbin = dir(fullfile( s.WORKBENCH, 'bin_*' ));
    if isscalar(wbbin)
        wbbin = fullfile( wbbin.folder, wbbin.name );
        setenv('PATH',sprintf('%s%s%s',wbbin,pathsep,getenv('PATH')));
        if isunix && system('hash wb_command')
            error('Workbench path is specified in osl.conf but it does not appear to be correct');
        end
    end

end