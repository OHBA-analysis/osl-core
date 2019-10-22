function s = read(fname)
%
% s = osl_conf.read( fname = $OSLCONF )
%
% Read config file.
%
% JH

    if nargin < 1, fname=getenv('OSLCONF'); end
    assert( ~isempty(fname), 'Missing path to config file.' );
    assert( osl_util.isfile(fname), 'File not found: %s', fname );

    % Parse the file
    f = fopen(fname,'r');
    l = fgetl(f);
    s = struct();
    while l ~= -1
        if ~isempty(strtrim(l))
            parts = strsplit(l,'=');
            name  = strtrim(parts{1});
            value = strjoin(parts(2:end),'=');
            s.(name) = strtrim(value);
        end
        l = fgetl(f);    
    end
    fclose(f);
    
    req = { 'FSLDIR', 'FSLBIN' ,'FSLLIB', 'SPMDIR', 'WORKBENCH' };
    assert( all(isfield(s,req)), 'Missing one of the following configuration fields:\n%s', strjoin(req,', ') );

end