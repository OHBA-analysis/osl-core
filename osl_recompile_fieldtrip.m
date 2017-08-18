function [replaced,notfound] = osl_recompile_fieldtrip( fieldtrip_dir )
    %
    % [replaced,notfound] = osl_recompile_fieldtrip( [fieldtrip_dir] )
    %
    % if fieldtrip_dir is omitted, it defaults to OSLDIR/spm12/external/fieldtrip.
    %
    %
    % There seem to be recurrent problems caused by Mex files in FieldTrip.
    % This script cleans up, compiles and replaces all occurrences of Mex files in FieldTrip.
    %
    % JH

    if ispc
        error('FT compilation via OSL is not supported on Windows');
    end

    OSLDIR = getenv('OSLDIR');

    current_dir = pwd;
    if nargin < 1
        fieldtrip_dir = fullfile( OSLDIR, 'spm12','external','fieldtrip' );
    end
    
    % delete all source Mex-files
    src_clean    = find_and_delete_mex(fullfile( fieldtrip_dir, 'src' ));
    fileio_clean = find_and_delete_mex(fullfile( fieldtrip_dir, 'fileio','@uint64' ));
    
    % extract function names to check for conflicts
    src_names    = cellfun( @basename, src_clean, 'UniformOutput', false );
    fileio_names = cellfun( @basename, fileio_clean, 'UniformOutput', false );
    assert( isempty(intersect(src_names,fileio_names)), 'Conflicts between Mex function names.' );
    
    % map all other Mex files
    mex_files = find_mex_files_recursive(fieldtrip_dir);
    mex_clean = clean_mex_files(mex_files);
    
    % recompile all source files
    cd(fieldtrip_dir);
    ft_compile_mex(true);
    cd(current_dir);
    
    % copy newly compiled version instead
    replaced = {};
    notfound = {};

    % Refresh the list of available src names
    src_names = union(src_names,cellfun( @basename, find_mex_files(fullfile(fieldtrip_dir,'src')), 'UniformOutput', false )); 

    for i = 1:numel(mex_clean) 
        
        file = mex_clean{i};
        name = basename(file);
        
        switch name
            case src_names
                delete([file '.mex*']);
                copyfile( fullfile(fieldtrip_dir,'src',[name '.' mexext]), fileparts(file) );
                replaced{end+1} = file;
            case fileio_names
                delete([file '.mex*']);
                copyfile( fullfile(fieldtrip_dir,'fileio','@uint64',[name '.' mexext]), fileparts(file) );
                replaced{end+1} = file;
            otherwise
                warning('Could not find a replacement for Mex-file: "%s"',mex_clean{i});
                notfound{end+1} = file;
        end
    end

    % Clean up src
    % Safer to leave disabled because one failure mode is if the original mex file is missing
    % in which case it won't be moved
    % delete(fullfile(fieldtrip_dir,'src',['*.',mexext]))
    
end

function b = basename(file)
    [~,b]=fileparts(file); 
end

function r = remext(file)
    [p,n]=fileparts(file);
    r = fullfile(p,n);
end

function mex_files = find_mex_files(folder)
    
    mex_files = dir(fullfile( folder, '*.mex*' ));
    mex_files = cellfun( @(x) fullfile(folder,x), {mex_files.name}, 'UniformOutput', false );
    
end

function mex_files = find_mex_files_recursive(folder)

    [~,mex_files] = system(['find "' folder '" -type f -name "*.mex*"']);
    mex_files = regexp(mex_files,'\n','split');
    mex_files = mex_files(~cellfun( @isempty, mex_files, 'UniformOutput', true ));

end

function mex_map = find_and_delete_mex(folder)

    mex_files = find_mex_files(folder);
    mex_map   = clean_mex_files(mex_files);
    
    if ~isempty(mex_files)
        for i = 1:numel(mex_files)
            delete(mex_files{i});
        end
    end

end

function mex_clean = clean_mex_files(mex_files)

    mex_clean = unique(cellfun( @remext, mex_files, 'UniformOutput', false ));
    
end
