function D = add_montage(D,W,name,labels,channels)
    % Add an online montage to an MEEG object (in memory)
    % 
    % INPUTS
    % - D (MEEG object)
    % - W (transformation matrix being applied to current montage (n_new_channels x n_old_channels)
    %     This matrix premultiplies the active montage and must therefore have the same number 
    % - name (optional string giving name of new montage)
    % - labels (optional cell array of names for each channel)
    % - channels (optional struct array that provides channel information)

    % First, check proposed tra matrix is OK in terms of orientation
    % That way, size(W,1) is guaranteed to be the new number of channels
    %
    % Romesh Abeysuriya 2017

    assert(size(W,2)==D.nchannels,sprintf('Transformation matrix is wrong size (is %dx%d, needs to be Nx%d)',size(W,1),size(W,2),D.nchannels));

    if nargin < 5 || isempty(channels) 
        channels = [];
    end
    
    if nargin < 4 || isempty(labels) 
        labels = arrayfun(@(x) sprintf('ROI%d',x),1:size(W,1),'UniformOutput',false);
    else
        assert(iscell(labels) && length(labels)==size(W,1),sprintf('Channel labels must be a cell array with %d elements',size(W,1)));
    end

    if nargin < 3 || isempty(name) 
        name = 'New montage';
    end

    currentMontage  = D.montage('getmontage');
    
    if isempty(currentMontage) % If there is currently no active montage - it's OK to set default channel information
        newMontage = struct; 
        newMontage.tra  = W;
        newMontage.labelorg = D.chanlabels;
        if ~isempty(channels) % If no channels provided AND no existing montage, SPM's default is fine
            newMontage.channels = channels;
        end
    else 
        % There are no original labels, and the tra matrix goes on top of the existing tra
        % The transformation is then from the currentMontage's labelorg to the new labels
        % But the channel types come from transforming the current montage channel types
        newMontage      = currentMontage;
        newMontage.tra  = W * currentMontage.tra;

        if ~isempty(channels)
            newMontage.channels = channels;
        else % If no channels provided AND no existing montage, SPM's default doesn't work because we need to read from the previous montage

            newMontage.channels = struct;

            for j = 1:size(W,1) % For each new channel

                old_channels = find(W(j,:)); % Old channels in the previous montage referred to by the new channel

                newMontage.channels(j).label = labels{j};

                % Bad if any old channels were bad
                newMontage.channels(j).bad = any([currentMontage.channels(old_channels).bad]); 

                % If only one previous chantype, use it - otherwise, call it 'Other'
                old_chantypes = unique({currentMontage.channels(old_channels).type});
                if numel(old_chantypes)==1
                    newMontage.channels(j).type = old_chantypes{1};
                else
                    % mixing different types of channels
                    newMontage.channels(j).type = 'Other';
                end

                % Same with units
                old_units = unique({currentMontage.channels(old_channels).units});
                if numel(old_units)==1
                    newMontage.channels(j).units = old_units{1};
                else
                    newMontage.channels(j).units = 'unknown';
                end
            end

        end

    end

    newMontage.name = name;
    newMontage.labelnew = labels;

    D = D.montage('add', newMontage);
    D = D.montage('switch', D.montage('getnumber'));
    