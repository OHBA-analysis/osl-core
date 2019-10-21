function save_figs(h,dirname,prefix)
    % Take in an array of figure handles
    % Save them to directory 'dirname'
    % Create directory if required
    % Automatically infer filename if required

    if nargin < 3 || isempty(prefix) 
        prefix = '';
    else
        prefix = [prefix '_'];
    end
    
    if ~isdir(dirname)
        mkdir(dirname);
    end

    for j = 1:length(h)
        fname = get(h(j),'tag');
        if isempty(fname)
            [~,tmp] = fileparts(tempname('.'));
            fname = ['untitled_' tmp(1:10)];
        end
        
        pfig(h(j),fullfile(dirname,[prefix fname '.png']));
    end
