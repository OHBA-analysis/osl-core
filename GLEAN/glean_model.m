function glean_model(GLEAN)
% Runs the model stage of GLEAN

if ~exist(GLEAN.model.model,'file')
    
    % Concatenate data:
    dataConcat = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
    subIndx    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
    for session = 1:numel(GLEAN.data)
        D = spm_eeg_load(GLEAN.subspace.data{session});
        dat = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
        dataConcat{session} = dat;
        subIndx{session} = session*ones(1,D.nsamples);

    end
    dataConcat = cell2mat(dataConcat);
    subIndx = cell2mat(subIndx); %#ok
    
    dataConcat = normalise(dataConcat,2); % TODO: maybe add an option for this
    
    switch lower(char(fieldnames(GLEAN.model.settings)))
        case 'hmm'
            hmmSettings = struct('K',GLEAN.model.settings.hmm.nstates,    ...
                                 'order',0,                               ...
                                 'Ninits',GLEAN.model.settings.hmm.nreps, ...
                                 'zeromean',0);
            hmm = glean_infer_hmm(dataConcat,hmmSettings); %#ok
            save(GLEAN.model.model,'hmm','subIndx')
        case 'ica'
            nICs = GLEAN.model.settings.ica.order;
            [ica.tICs,ica.SM,~] = fastica(dataConcat,     ...
                                          'g','tanh',     ...
                                          'lastEig',nICs, ...
                                          'numOfIC',nICs, ...
                                          'approach','symm');
            save(GLEAN.model.model,'ica','subIndx')
    end
    
end

end