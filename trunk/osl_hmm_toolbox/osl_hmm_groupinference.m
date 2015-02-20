function [HMMresults,statemaps,epoched_statepath_sub] = osl_hmm_groupinference(data_files,hmmdir,todo,options)
% Runs a group HMM analysis on amplitude envelopes of source reconstructed
% MEG data using the SPM12 beamformer OR any arbitrary data saved in .mat
% data files.
%
% osl_hmm_groupinference(data_files,hmmdir,todo,env_settings,hmm_settings)
%
% INPUTS:
%
% data_files - cell array containing list of SPM12 source recon results 
%              containing continuous data OR .mat files containing data 
%              [observations x timepoints]
%  
% hmmdir     - directory within which to save HMM results
%
% todo       - binary flags for each stage of the analysis [0/1]:
%              .envelope - computes amplitude envelopes from source
%                        reconstructed data
%              .concat   - concatenates data and performs PCA dimensionality
%                        reduction
%              .hmm      - infers the HMM
%              .output   - output the state maps
%
%
% options    - struct with settings for envelope, concatenation and hmm stages
%              with the following fields:
%              .envelope - settings for envelope computation with fields:
%                          .windowsize - window size for moving average filter [s]
%                            (default 0.1)
%                          .multiband  - multiband frequencies {[Hz Hz],[Hz Hz]}
%                          .filenames  - filenames to save/load envelopes from
%                            (default <hmmdir>/envelopes/<BFfiles{subnum}>.mat)
%             
%              .concat   - settings for temporal concatenation with fields:
%                          .log        - apply log transform
%                          .norm_subs  - apply subject variance normalisation [0/1]
%                            (default 1)
%                          .pcadim     - dimensionality to use 
%                            (default 40)                       
%                          .whiten     - apply whitening [0/1]                      
%                            (default 1)
%                          .filename   - filename to save/load concatenated data from
%                            (default <hmmdir>/env_concat.mat)
%
%              .hmm      - settings for HMM inference with fields:                
%                          .nstates    - number of states to infer 
%                            (default 8)
%                          .nreps      - number of repeat inferences to run 
%                            (default 5)
%                          .filename   - filename to save/load HMM results from
%                            (default <hmmdir>/HMM.mat)
%      
%              .output   - settings for the HMM output with fields
%                          .method     - maps/matrix to output:
%                                        'pcorr' - partial correlation maps
%                                        'var'   - variance maps
%                                        'cov'   - covariance matrices
%                            (default 'pcorr')
%                          .filename   - filename to save/load HMM results from
%                            (default <hmmdir>/HMM<_method.*>)
%                       
%
% OUTPUTS:
%
% HMMresults - filename of the inferred HMM results
%
% statemaps  - filename of the HMM state maps
%
%
% AB 2014

% global OSLDIR

HMMresults = [];
statemaps  = [];

data_files = cellstr(data_files);

data_fnames = cell(size(data_files));
for f = 1:numel(data_files)
  [~,data_fnames{f},~] = fileparts(data_files{f});
end

if ~isdir(hmmdir)
  mkdir(hmmdir); 
end


% Set up directories:
if isfield(options.envelope,'filenames') && ~isempty(options.envelope.filenames)
    filenames.envelope = cell(size(options.envelope.filenames));
    [pathstr,filestr] = cellfun(@fileparts,options.envelope.filenames,'uniformoutput',0);
    for f = 1:length(options.envelope.filenames)
        if strcmp(pathstr(f),'')
            pathstr{f} = fullfile(hmmdir,'envelopes');
            if ~isdir(pathstr{f})
                mkdir(pathstr{f})
            end
        end
        filenames.envelope{f} = fullfile(pathstr{f},filestr{f});
    end
else
    pathstr = fullfile(hmmdir,'envelopes');
    if ~isdir(pathstr)
        mkdir(pathstr)
    end
    filenames.envelope = fullfile(pathstr,data_fnames);
end
filenames.envelope = strcat(filenames.envelope, '.mat');

if isfield(options.concat,'filename') && ~isempty(options.concat.filename)
    [pathstr,filestr] = fileparts(options.concat.filename);
        if strcmp(pathstr,'')
            pathstr = hmmdir;
        end
        filenames.concat = fullfile(pathstr,filestr);
else
    filenames.concat = fullfile(hmmdir,'env_concat');
end
filenames.concat = [filenames.concat '.mat'];

if isfield(options.hmm,'filename') && ~isempty(options.hmm.filename)
    [pathstr,filestr] = fileparts(options.hmm.filename);
        if strcmp(pathstr,'')
            pathstr = hmmdir;
        end
        filenames.hmm = fullfile(pathstr,filestr);
else
    filenames.hmm = fullfile(hmmdir,'hmm');
end
filenames.hmm = [filenames.hmm '.mat'];

if isfield(options.output,'filename') && ~isempty(options.output.filename)
    [pathstr,filestr] = fileparts(options.output.filename);
        if strcmp(pathstr,'')
            pathstr = hmmdir;
        end
        filenames.output = fullfile(pathstr,filestr);
else
    filenames.output = fullfile(hmmdir,'hmm');
end



% Default todo settings
try todo.envelope = todo.envelope; catch, todo.envelope = 1; end
try todo.concat   = todo.concat;   catch, todo.concat   = 1; end
try todo.hmm      = todo.hmm;      catch, todo.hmm      = 1; end

% Default envelope settings
try windowsize = options.envelope.windowsize; catch, windowsize = 0.1;  end
try multiband  = options.envelope.multiband;  catch, multiband  = [];   end
try envelope_do  = options.envelope.do;  catch, envelope_do  = 1;   end

% Default concatenation settings
try logtrans  = options.concat.log;        catch, logtrans   = 0;  end
try norm_subs = options.concat.norm_subs;  catch, norm_subs  = 1;  end
try pcadim    = options.concat.pcadim;     catch, pcadim     = 40; end
try whiten    = options.concat.whiten;     catch, whiten     = 1;  end

% Default HMM settings
try nstates = options.hmm.nstates; catch, nstates = 8; end
try nreps   = options.hmm.nreps;   catch, nreps   = 5; end
    
% Default output settings
try output_method = options.output.method; catch, output_method = 'pcorr'; end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    C O M P U T E   E N V E L O P E S                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.envelope
  
  for subnum = 1:length(data_files)
    
    [pathstr,filestr] = fileparts(data_files{subnum});
  
    disp(['Computing envelope data for ' filestr]);
   
    % ---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---
    if ~isempty(multiband)
        for f = 1:numel(multiband)
            D = spm_eeg_load(data_files{subnum});
            current_montage = montage(D,'getindex');
            D = montage(D,'switch',0);
            spm_eeg_filter(struct('D',D,'band','bandpass','freq',multiband{f}));
            D = spm_eeg_load(prefix(data_files{subnum},'f'));
            D = montage(D,'switch',current_montage); D.save;
            
            S = [];
            S.D = fullfile(D.path,D.fname);
            S.winsize = windowsize;
            
            D = osl_hilbenv(S);
            move(D,strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            D = spm_eeg_load(prefix(data_files{subnum},'f'));
            D.delete;
            disp(['Saving envelope data for ' filestr ' to ' filenames.envelope{subnum} ' band ' num2str(f)])
            clear D
        end
    else % Single frequency band
        D = spm_eeg_load(data_files{subnum});
        S = [];
        S.D = fullfile(D.path,D.fname);
        %S.winsize = fix(windowsize*D.fsample);
        S.winsize = windowsize;
        D = montage(D,'switch',2); D.save;
        
        S2=[];
        S2.D=S.D;
        S2.outfile=prefix(D.fnamedat,'raw');
        Draw=spm_eeg_copy(S2);

        if envelope_do
            Denv = osl_hilbenv(S);          
            move(Denv,filenames.envelope{subnum});  
            disp(['Saving envelope data for ' filestr ' to ' filenames.envelope{subnum}])
        end;
        
        [pth fname ext]=fileparts(filenames.envelope{subnum});
        newrawname=[pth '/' fname '_raw' ext];
        move(Draw,newrawname);  
        disp(['Saving raw data for ' filestr ' to ' newrawname])

        clear D Denv Draw
    end
    % ---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---
    
    
  end
  
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P A R E   H M M   D A T A                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.concat || (todo.infer && ~exist(filenames.concat,'file'))
    
    for f = 1:max(numel(multiband),1)
        
        % Load subjects and concatenate:
        dat_concat = [];
        subj_inds=[];
        C = 0;
        
        for subnum = 1:length(data_files)
            
            if ~isempty(multiband)
                D = spm_eeg_load(strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            else
                [pth fname ext]=fileparts(filenames.envelope{subnum});
                newrawname=[pth '/' fname '_raw' ext];      
                Draw = spm_eeg_load(newrawname);
                
                if envelope_do
                    Denv = spm_eeg_load(filenames.envelope{subnum});
                    Ddat = Denv;
                    Dcov = Denv;
                else
                    Ddat = Draw;
                    Dcov = Draw;                    
                end;
                clear Denv D Draw;
                
            end
          
            % use all conditions
            for condnum = 1:length(Ddat.condlist),
                
                trials = Ddat.indtrial(Ddat.condlist{condnum},'good'); 

                for trl=1:length(trials),                            
                                       
                    tbad = all(badsamples(Ddat,':',':',trials(trl)));
                    samples2use = find(~tbad);
                    dat = Ddat(:,samples2use,trials(trl)); %#ok - SPM doesn't like logical indexing

                    tbad = all(badsamples(Dcov,':',':',trials(trl)));
                    samples2use = find(~tbad);
                    cov_dat = Dcov(:,samples2use,trials(trl)); %#ok - SPM doesn't like logical indexing

                    if logtrans && envelope_do
                        dat = log10(dat);
                    end
                    
                    if norm_subs
                       
                        if envelope_do
                            dat = demean(dat,2)./std(dat(:));
                            %cov_dat = demean(cov_dat,2)./std(cov_dat(:));
                            
                            %dat = normalise(dat,2);
                            cov_dat = normalise(cov_dat,2);
                        else
                            %dat = demean(dat,2)./std(dat(:));
                            %cov_dat = demean(cov_dat,2)./std(cov_dat(:));
                    
                            dat = normalise(dat,2);
                            cov_dat = normalise(cov_dat,2);
                            
                        end;
                    end

                    dat_concat = [dat_concat,dat];

                    subj_inds = [subj_inds,subnum*ones(1,size(cov_dat,2))];            
                    %cond_inds = [cond_inds,condnum*ones(1,size(dat,2))];

                    C = C + cov_dat*permute(cov_dat,[2,1]); 
                end;
            end;
        end
        C = C ./ (length(subj_inds)-1);
        clear cov_dat dat subj_inds;
        
        % PCA + whitening - below is equivalent to fastica whitening code but much faster
        [allsvd,MixingMatrix] = eigdec(C,pcadim);
        if whiten
            MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
        end
        hmmdata =  (MixingMatrix * dat_concat)';
        fsample = Ddat.fsample;
        
        if ~isempty(multiband)
            savestr = strrep(filenames.concat,'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']);
            disp(['Saving concatenated envelope data to ' savestr])
            save(savestr,'hmmdata','MixingMatrix','fsample')
        else
            disp(['Saving concatenated envelope data to ' filenames.concat])
            save(filenames.concat,'hmmdata','MixingMatrix','fsample')
        end
        
        % save mixing matrix to nii
        [pth nm ext]=fileparts(filenames.concat);
        pc_maps=[pth '/' nm '_pc_maps'];
        pc_maps = nii_quicksave(MixingMatrix',pc_maps,getmasksize(Dcov.nchannels),2);
        disp(['PCA maps saved to ' pc_maps]);   
        
        [pth nm ext]=fileparts(filenames.concat);
        mean_pc_maps=[pth '/' nm '_mean_pc_maps'];
        mean_pc_maps = nii_quicksave(mean(abs(MixingMatrix),1)',mean_pc_maps,getmasksize(Dcov.nchannels),2);
        disp(['mean of PCA maps saved to ' mean_pc_maps]);  
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            I N F E R   H M M                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.infer
    
    if ~isempty(multiband)
        hmmdata = [];
        freq_inds = [];
        MixingMatrix = cell(numel(multiband),1);
        for f = 1:numel(multiband)
            tmp = load(strrep(filenames.concat,'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            hmmdata = [hmmdata tmp.hmmdata];
            MixingMatrix{f} = tmp.MixingMatrix;
        end

        fsample = tmp.fsample;
    else
        load(filenames.concat)
    end
    
    %%%%%%%%%%%%%%%%
    % Sort this OUT!
    if 0
        %tilde='/Users/woolrich';
        %addpath(genpath([tilde '/homedir/vols_data/ctf_selfpaced_fingertap']));
        global OSLDIR
        rmpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox')));
        addpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

        if ~envelope_do
            hmm = ABhmm_infer(hmmdata,nstates,nreps,'constrain_mean');
        else
            hmm = ABhmm_infer(hmmdata,nstates,nreps);
        end;
        hmm.statepath = ABhmm_statepath(hmm);
        addpath(genpath(OSLDIR));

    else
        global OSLDIR
        rmpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
        addpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox')));
        
        hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',~envelope_do));
        addpath(genpath(OSLDIR));
    end;
    %%%%%%%%%%%%%%%%
    
    hmm.MixingMatrix = MixingMatrix;
    hmm.fsample = fsample;

    
    % Save results
    disp(['Saving inferred HMM to ' filenames.hmm])    
    save(filenames.hmm,'hmm');
    
    HMMresults = filenames.hmm;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       O U T P U T   R E S U L T S                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.output
    
    for f = 1:max(numel(multiband),1)
        
        if ~isempty(multiband)
            statemaps = [filenames.output,'_',output_method,'_',strrep(num2str(multiband{f}),'  ','_') 'Hz'];
        else
            statemaps = [filenames.output,'_',output_method];
        end
        
        if ~exist('hmm','var')
            HMMresults = filenames.hmm;
            load(HMMresults)
        end
        
        switch output_method
            case {'pcorr','var'}
                % Compute partial correlation/variance maps
                if ~isempty(multiband)
                    stat = zeros(size(hmm.MixingMatrix{1},2),hmm.K);
                else
                    stat = zeros(size(hmm.MixingMatrix,2),hmm.K);
                end
                
                subjstart_index=1;
                
                %%%% find out number of conditions
                subnum=1;
                if ~isempty(multiband)
                    D = spm_eeg_load(strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
                else
                    if envelope_do
                        D = spm_eeg_load(filenames.envelope{subnum});
                    else
                        [pth fname ext]=fileparts(filenames.envelope{subnum});
                        newrawname=[pth '/' fname '_raw' ext];      
                        D = spm_eeg_load(newrawname);
                    end;             
                end
                num_conds=length(D.condlist);
                clear D;
                %%%%%
                
                epoched_statepath_sub=cell(length(data_files),num_conds);
                
                for subnum = 1:length(data_files)
                    
                    %try
                        if ~isempty(multiband)
                            D = spm_eeg_load(strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
                        else
                            if envelope_do
                                D = spm_eeg_load(filenames.envelope{subnum});
                            else
                                [pth fname ext]=fileparts(filenames.envelope{subnum});
                                newrawname=[pth '/' fname '_raw' ext];      
                                D = spm_eeg_load(newrawname);                
                            end;
                              
                        end
                        disp(['Computing ' output_method ' maps for ' D.fname]);
                        
                        dat_concat_sub=[];                                                
                        
                        % use all conditions
                        trialstart_index=1;
                        for condnum = 1:num_conds,
                            trials = D.indtrial(D.condlist{condnum},'good'); 
                            
                            epoched_statepath_sub{subnum, condnum}=zeros(1,D.nsamples,length(trials));
                                                      
                            for trl=1:length(trials),                            
                                tbad = all(badsamples(D,':',':',trials(trl)));
                                
                                samples2use = find(~tbad);
                                dat = D(:,samples2use,trials(trl)); %#ok - SPM doesn't like logical indexing                                                                
                                
                                if logtrans
                                    dat = log10(dat);
                                end
                                if norm_subs
                                    if envelope_do
                                        %dat = normalise(dat,2);                                        
                                        dat = demean(dat,2)./std(dat(:));
                                    else
                                        dat = normalise(dat,2);   
                                        %dat = demean(dat,2)./std(dat(:));
                                    end;
                                end

                                dat_concat_sub = [dat_concat_sub,dat];                                                                    
                                
                                epoched_statepath_sub{subnum, condnum}(1,samples2use,trl)=hmm.statepath(trialstart_index:trialstart_index+size(dat,2)-1);
                                trialstart_index=trialstart_index+size(dat,2);
                                                    
                            end;
                        end;                    
                    
                        hmm_sub = hmm; 
                        from=subjstart_index;
                        to=from+(trialstart_index-1)-1;
                        hmm_sub.statepath = hmm.statepath(from:to);                    
                        subjstart_index=to+1;
                        
                        stat = stat + osl_hmm_statemaps(hmm_sub,dat_concat_sub,~envelope_do,output_method);
                    
                    %catch
                    %    error([output_method ' currently only supported for OAT results'])
                    %end
                                
                end
                
                stat = stat ./ length(data_files);
                disp(['Saving state spatial maps to ' statemaps])
                statemaps = nii_quicksave(stat,statemaps,getmasksize(D.nchannels),2);
               
                % resave updated hmm
                hmm.statemaps=statemaps;
                hmm.epoched_statepath_sub=epoched_statepath_sub;
                disp(['Saving updated hmm to ' filenames.hmm]);
                save(filenames.hmm,'hmm');        
                
            case 'cov'
                cov_mats = osl_hmm_statemaps(hmm,[],[],'cov');
                statemaps = [statemaps '.mat'];
                disp(['Saving state covariance matrices to ' statemaps])
                save(statemaps,'cov_mats')
                
            otherwise
                warning([output_method ' is not a supported output method']);
        end
        
    end
    
end


