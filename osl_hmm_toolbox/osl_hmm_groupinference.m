function [HMMresults,statemaps] = osl_hmm_groupinference(data_files,hmmdir,todo,options)
% Runs a group HMM analysis on amplitude envelopes of source reconstructed
% MEG data using the SPM12 beamformer OR any arbitrary data saved in .mat
% data files.
%
% osl_hmm_groupinference(data_files,hmmdir,todo,options)
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

OSLDIR = getenv('OSLDIR');


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

% Default prepare settings
try envelope_do = options.prepare.envelope;   catch, envelope_do = 1;   end

% Default envelope settings
try windowsize = options.envelope.windowsize; catch, windowsize = 0.1;  end
try multiband  = options.envelope.multiband;  catch, multiband  = [];   end

% Default concatenation settings
try logtrans  = options.concat.log;        catch, logtrans   = 0;  end
try norm_subs = options.concat.norm_subs;  catch, norm_subs  = 1;  end
try pcadim    = options.concat.pcadim;     catch, pcadim     = 40; end
try whiten    = options.concat.whiten;     catch, whiten     = 1;  end

% Default HMM settings
try nstates = options.hmm.nstates; catch, nstates = 8; end
try nreps   = options.hmm.nreps;   catch, nreps   = 5; end
try use_old_tbx = options.hmm.use_old_hmm_tbx;  catch, use_old_tbx = 0; end
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
        S.winsize = windowsize;
        
        D = osl_hilbenv(S);
        move(D,filenames.envelope{subnum});  
        disp(['Saving envelope data for ' filestr ' to ' filenames.envelope{subnum}])
        clear D
    end
    % ---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---
    
    
  end
  
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P A R E   H M M   D A T A                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.concat || (todo.infer && ~osl_util.isfile(filenames.concat))
    
    for f = 1:max(numel(multiband),1)
        
        % Load subjects and concatenate:
        env_concat = [];
        subj_inds = [];
        C = 0;
        
        for subnum = 1:length(data_files)
            
            if ~isempty(multiband)
                D = spm_eeg_load(strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            else
                D = spm_eeg_load(filenames.envelope{subnum});
            end
            tbad = ~good_samples(D);
            samples2use = find(~tbad);
            env = D(:,samples2use); %#ok - SPM doesn't like logical indexing
            if logtrans
                env = log10(env);
            end
            if norm_subs
                env = demean(env,2)./std(env(:));
                %env = normalise(env,2);
            end
            
            env_concat = [env_concat,env];
            subj_inds = [subj_inds,subnum*ones(1,size(env,2))];
            C = C + env*permute(env,[2,1]);
        end
        C = C ./ (length(subj_inds)-1);
        clear env
        
        
        % PCA + whitening - below is equivalent to fastica whitening code but much faster
        [allsvd,MixingMatrix] = eigdec(C,pcadim);
        if whiten
            MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
        end
        hmmdata =  (MixingMatrix * env_concat)';
        fsample = D.fsample;
        
        if ~isempty(multiband)
            savestr = strrep(filenames.concat,'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']);
            disp(['Saving concatenated envelope data to ' savestr])
            save(savestr,'hmmdata','MixingMatrix','fsample','subj_inds')
        else
            disp(['Saving concatenated envelope data to ' filenames.concat])
            save(filenames.concat,'hmmdata','MixingMatrix','fsample','subj_inds')
        end
        
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
        subj_inds = tmp.subj_inds;
        fsample = tmp.fsample;
    else
        load(filenames.concat)
    end
    
    
    %hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',0));
    
    
    % Switch between Iead's & Diego's HMM toolboxes
    if use_old_tbx
        rmpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
        addpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

        if envelope_do 
            hmm = ABhmm_infer(hmmdata,nstates,nreps);
        else
            hmm = ABhmm_infer(hmmdata,nstates,nreps,'constrain_mean');
        end;
        addpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
        rmpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

    else

        rmpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
        addpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));

        hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',~envelope_do));
        %hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',false));
        addpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
        rmpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));

    end
    
    
    hmm.MixingMatrix = MixingMatrix;
    hmm.fsample = fsample;
    
    % Save results
    disp(['Saving inferred HMM to ' filenames.hmm])
    save(filenames.hmm,'hmm','subj_inds')
    
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
                for subnum = 1:length(data_files)
                    
                    try
                        if ~isempty(multiband)
                            D = spm_eeg_load(strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
                        else
                            D = spm_eeg_load(filenames.envelope{subnum});
                        end
                        disp(['Computing ' output_method ' maps for ' D.fname]);
                        
                        tbad = ~good_samples(D);
                        samples2use = find(~tbad);
                        env = D(:,samples2use); %#ok - SPM doesn't like logical indexing
                    catch
                        error([output_method ' currently only supported for OAT results'])
                    end
                    
                    hmm_sub = hmm; hmm_sub.statepath = hmm.statepath(subj_inds==subnum);
                    stat = stat + osl_hmm_statemaps(hmm_sub,env,0,output_method);
                end
                stat = stat ./ length(unique(subj_inds));
                disp(['Saving state spatial maps to ' statemaps])
                statemaps = nii.quicksave(stat,statemaps,getmasksize(D.nchannels),2);
                
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

