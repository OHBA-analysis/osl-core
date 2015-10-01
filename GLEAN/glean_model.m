function glean_model(GLEAN)
% Runs the model stage of GLEAN

if ~exist(GLEAN.model.model,'file')
    
    % Concatenate data:
    dataConcat = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
    subIndx    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
    for session = 1:numel(GLEAN.data)
        D = spm_eeg_load(GLEAN.subspace.data{session});
        %if isfield(GLEAN.settings.envelope,'freqbands') && numel(GLEAN.settings.envelope.freqbands) > 1
        %    % rearrange channels x frequencies as [c1f1, ... ,c1fn, c2f1, ...]
            dat = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
            
            %dataConcat{session} = normalise(dat,2);
            dataConcat{session} = dat;
        %else
        %    dataConcat{session} = D(:,:,:);
        %end
        
        % Commented code above should be deleted once 4D MEEG object is
        % used consistently
        
        subIndx{session} = session*ones(1,D.nsamples);

    end
    dataConcat = cell2mat(dataConcat);
    subIndx = cell2mat(subIndx); %#ok
    
    dataConcat = normalise(dataConcat,2); % MOVE THIS SOMEWHERE BETTER!
    
    switch lower(char(fieldnames(GLEAN.model.settings)))
        case 'hmm'
            hmm = osl_hmm_infer(dataConcat,struct('K',GLEAN.model.settings.hmm.nstates,'order',0,'Ninits',GLEAN.model.settings.hmm.nreps,'zeromean',0)); %#ok
            save(GLEAN.model.model,'hmm','subIndx')
        case 'ica'
            nICs = GLEAN.model.settings.ica.order;
            [ica.tICs,ica.SM,~] = fastica(dataConcat,'g','tanh','lastEig',nICs,'numOfIC',nICs,'approach','symm');
            save(GLEAN.model.model,'ica','subIndx')
    end
    
end

end