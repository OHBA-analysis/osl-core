function varargout = has_montage(D,pattern,case_sensitive)
    % Check whether an MEEG has montages with names satisfying a search
    %
    % [n_present,montage_number] = has_montage(D,pattern,case_sensitive)
    %
    % INPUTS
    % - D - An MEEG object
    % - pattern - A search string (default='')
    % - case_sensitive - If true, case must match as well (default=false)
    %
    % If you run 'has_montage(D)' with no other arguments, a summary of the
    % available montages will be printed. An asterisk denotes the active montage.
    %
    % OUTPUTS
    % - n_present - Number of montages containing the pattern in the montage name
    % - montage_number - List of length n_present specifying the montage index of each match
    %
    % TYPICAL USAGE
    % 
    % - If you know that your analysis pipeline makes a montage with a particular
    %   name, you can check if that name is present to see if the pipeline has been run
    % - If your function requires a particular montage to be selected, you can find
    %   the required montage number using this function, and then switch to it 
    %
    % EXAMPLE
    %
    % >> has_montage(D)
    % 0 - none (276 channels)
    % 1 - without weights normalisation, class 1 (3559 channels)
    % 2 - with weights normalisation, class 1 (3559 channels)
    % *3 - Parcellated with weights normalisation, class 1 (38 channels)
    % >> [n,idx] = has_montage(D,'parcellated')
    %   n =
    %        1
    %   idx =
    %        3
    % >> D.montage('getmontage',idx).name
    %   ans =
    %       'Parcellated with weights normalisation, class 1'

    if nargin < 3 || isempty(case_sensitive) 
        case_sensitive = false;
    end

    if nargin < 2 || isempty(pattern) 
        pattern = [];
    end
    
    n_present = 0;
    montage_number = [];
    n_montages = D.montage('getnumber');
    active_montage = D.montage('getindex');

    if nargout > 0
        varargout{1} = n_present;
        varargout{2} = montage_number;
    end
    
    if isempty(pattern)
        % If no search pattern, print the montages and return
        for j = 0:n_montages
            D = D.montage('switch',j);
            if active_montage == j
                fprintf('*');
            end
            fprintf('%d - %s (%d channels)\n',j,D.montage('getname'),D.nchannels);
        end
        return
    end
    
    if n_montages == 0 % No online montages
        return
    end

    for j = 1:n_montages
        m = D.montage('getmontage',j);

        if case_sensitive
            f = regexp(m.name,pattern,'match');
        else
            f = regexpi(m.name,pattern,'match');
        end

        if ~isempty(f)
            n_present = n_present + 1;
            montage_number(end+1) = j;
        end
    end

    if nargout > 0
        varargout{1} = n_present;
        varargout{2} = montage_number;
    end