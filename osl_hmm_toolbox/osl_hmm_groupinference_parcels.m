function [HMMresults,statemaps,epoched_statepath_sub] = osl_hmm_groupinference_parcels(data_files,hmmdir,todo,options)
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
%              .prepare - preps data, including computes amplitude prepares from source
%                        reconstructed data
%              .concat   - concatenates data and performs PCA dimensionality
%                        reduction
%              .hmm      - infers the HMM
%              .output   - output the state maps
%
%
% options    - struct with settings for preparing data and hmm stages
%              with the following fields:
%              .prepare - settings for prepare computation with fields:
%                          .envelope   - compute Hilbert envelope for data
%                            (default 1)
%                          .windowsize - window size for moving average filter [s]
%                            only carried out if data is enveloped
%                            (default 0.1)
%                          .filenames  - filenames to save/load prepares from
%                            (default <hmmdir>/prepare/<BFfiles{subnum}>.mat)
%                          .parcellation.file - 3D niftii file with same
%                            gridstep as data_files containing parcellation
%                            labels. If not provided then no parcellation is
%                            carried out (default)
%                          .parcellation.method - passed to ROInets.get_corrected_node_tcs
%                          .parcellation.protocol - passed to ROInets.get_corrected_node_tcs
%                          .log - apply log transform to envelopes
%
%              .concat   - settings for temporal concatenation with fields:                      
%                          .pcadim        - dimensionality to use (default 40)             
%                          .normalisation - variance normalisation method:
%                            'none'     [data = demean(data,2)]
%                            'global'   [data = demean(data,2)/std(data(:)]
%                            'voxelwise'[data = normalise(data,2)]
%                            (default 'global')
%                          .whiten        - apply whitening [0/1] (default 1)
%                          .filename      - filename to save/load concatenated data from
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
%                                        
%                            (default 'pcorr')
%                          .use_parcel_weights - uses parcel weights rather
%                          than binary parcels [0/1] (default 1)
%                          
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
% AB & MWW 2014 

global OSLDIR

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
if isfield(options.prepare,'filenames') && ~isempty(options.prepare.filenames)
    filenames.prepare = cell(size(options.prepare.filenames));
    [pathstr,filestr] = cellfun(@fileparts,options.prepare.filenames,'uniformoutput',0);
    for f = 1:length(options.prepare.filenames)
        if strcmp(pathstr(f),'')
            pathstr{f} = fullfile(hmmdir,'prepare');
            if ~isdir(pathstr{f})
                mkdir(pathstr{f})
            end
        end
        filenames.prepare{f} = fullfile(pathstr{f},filestr{f});
    end
else
    pathstr = fullfile(hmmdir,'prepare');
    if ~isdir(pathstr)
        mkdir(pathstr)
    end
    filenames.prepare = fullfile(pathstr,data_fnames);
end
filenames.prepare = strcat(filenames.prepare, '.mat');

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
try todo.prepare = todo.prepare; catch, todo.prepare = 1; end
try todo.concat  = todo.concat;  catch, todo.concat  = 1; end
try todo.hmm     = todo.hmm;     catch, todo.hmm     = 1; end

% Default prepare settings
try envelope_do = options.prepare.envelope;   catch, envelope_do = 1;   end
try logtrans    = options.prepare.log;        catch, logtrans    = 0;   end
try windowsize  = options.prepare.windowsize; catch, windowsize  = 0.1; end

parcellation = [];
try parcellation.file     = options.prepare.parcellation.file;     catch, parcellation.file     = [];          end
try parcellation.protocol = options.prepare.parcellation.protocol; catch, parcellation.protocol = 'symmetric'; end
try parcellation.method   = options.prepare.parcellation.method;   catch, parcellation.method   = 'PCA';       end
use_parcels               = ~isempty(parcellation.file);

% Default concatenation settings
try pcadim        = options.concat.pcadim;        catch, pcadim         = 40;       end
try whiten        = options.concat.whiten;        catch, whiten         = 1;        end
try normalisation = options.concat.normalisation; catch, normalisation  = 'global'; end

% Default HMM settings
try nstates = options.hmm.nstates; catch, nstates = 8; end
try nreps   = options.hmm.nreps;   catch, nreps   = 5; end
    
% Default output settings
try output_method       = options.output.method;                catch, output_method = 'pcorr'; end
try use_parcel_weights  = options.output.use_parcel_weights;    catch, use_parcel_weights = 0;  end
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P   D A T A
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.prepare
    
    for subnum = 1:length(data_files)
        
        % Compute envelopes for all voxels
        if envelope_do
            S         = [];
            S.D       = data_files{subnum};
            S.winsize = windowsize;
            D = osl_hilbenv(S);
        else
            D = spm_eeg_load(data_files{subnum});
        end
        
        % Copy MEEG object to HMM directory
        move(D,filenames.prepare{subnum}); 
        disp(['Saving envelope data for ' filestr ' to ' filenames.prepare{subnum}])
        clear D
        
        
        % Apply parcellation
        if use_parcels
            
            S                   = [];
            S.D                 = data_files{subnum};
            S.parcellation      = parcellation.file;
            S.orthogonalisation = parcellation.protocol;
            S.method            = parcellation.method;
            S.prefix            = 'p'; 
            [D,parcelWeights,parcelAssignments] = osl_apply_parcellation(S);
            
            % Save parcellation results:
            pathstr = fileparts([filenames.prepare{1}]);
            fname = fullfile(pathstr,'parcellation');
            save(fname,'parcelAssignments','parcelWeights');
            
            % Compute envelopes for all parcels
            if envelope_do
                S         = [];
                S.D       = fullfile(D.path,D.fname);
                S.winsize = windowsize;
                D = osl_spmfun(@osl_hilbenv,S); % keep same filename
            end
            
            D = move(D,[fileparts(filenames.prepare{1}),filesep]);
        
        end

    end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    C O N C A T   H M M   D A T A                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.concat || (todo.infer && ~exist(filenames.concat,'file'))
            
    % Load subjects and concatenate:
    dat_concat  = [];
    subj_inds   = [];
    C           = 0;
    Csamples    = 0;
    
    for subnum = 1:length(data_files)

        if use_parcels
            D = spm_eeg_load(prefix(filenames.prepare{subnum},'p'));
        else
            D = spm_eeg_load(filenames.prepare{subnum});
        end
        
        [data,cov_data] = prepare_data(D,normalisation,logtrans);
      
        % compute covariance
        C = C + cov_data * permute(cov_data,[2,1]);
        
        % keep track of number of samples
        Csamples = Csamples + size(cov_data,2);
        
        dat_concat = [dat_concat, data];
        subj_inds  = [subj_inds, subnum*ones(1,size(data,2))];

    end
    
    C = C / (Csamples-1);
    
    clear data cov_data;

    % Apply PCA dimensionality reduction
    if pcadim > 0
        disp(['Keeping top ' num2str(pcadim) ' PCs']);
        pcadim = min(pcadim,size(dat_concat,1));   
        [allsvd,MixingMatrix] = eigdec(C,pcadim);
    else
        MixingMatrix = eye(size(dat_concat,1));   
        allsvd = diag(C);
    end
    
    if whiten
        MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
    else
        MixingMatrix = MixingMatrix';
    end
    
    hmmdata = (MixingMatrix * dat_concat)';
                 
    fsample = D.fsample;
   
    disp(['Saving group concatenated data to ' filenames.concat])
    save(filenames.concat,'hmmdata','MixingMatrix','fsample','subj_inds');

    if use_parcels
        % Load parcellation results:
        pathstr = fileparts([filenames.prepare{1}]);
        fname   = fullfile(pathstr,'parcellation');
        load(fname); % Loads parcelAssignments and parcelWeights
        masksize = getmasksize(size(parcelAssignments,1));
    else
        masksize=getmasksize(size(dat_concat,1));
    end
    
    
    
    % Save PCA maps
    [pathstr,filestr] = fileparts(filenames.concat);
    savefile_pc_maps        = [pathstr '/' filestr '_pc_maps'];
    savefile_mean_pc_maps   = [pathstr '/' filestr '_mean_pc_maps'];
    savefile_std_maps    = [pathstr '/' filestr '_hmmdata_std_maps'];
    
    if use_parcels
        ROInets.nii_parcel_quicksave(MixingMatrix',parcelAssignments,savefile_pc_maps,masksize,masksize,'nearestneighbour');
        ROInets.nii_parcel_quicksave(mean(abs(MixingMatrix),1)',parcelAssignments,savefile_mean_pc_maps,masksize,2,'nearestneighbour');
        ROInets.nii_parcel_quicksave(std(pinv(MixingMatrix)*hmmdata',[],2),parcelAssignments,savefile_std_maps,masksize,2,'nearestneighbour');
    else
        nii_quicksave(MixingMatrix',savefile_pc_maps,masksize,2); 
        nii_quicksave(mean(abs(MixingMatrix),1)',savefile_mean_pc_maps,masksize,masksize);
        nii_quicksave(std(pinv(MixingMatrix)*hmmdata',[],2),savefile_std_maps,masksize,masksize);
    end

    disp(['PCA maps saved to ' savefile_pc_maps]);   
    disp(['mean of PCA maps saved to ' savefile_mean_pc_maps]);  
    disp(['HMM data sd maps saved to ' savefile_std_maps]);     

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            I N F E R   H M M                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.infer
    
    load(filenames.concat);
    D = spm_eeg_load(data_files{1}); 

    % Switch between Iead's & Diego's HMM toolboxes
    if options.hmm.use_old_hmm_tbx
        rmpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
        addpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

        if ~envelope_do
            hmm = ABhmm_infer(hmmdata,nstates,nreps,'constrain_mean');
        else
            hmm = ABhmm_infer(hmmdata,nstates,nreps);
        end;
        hmm.statepath = ABhmm_statepath(hmm);
        addpath(genpath(OSLDIR));

    else
        rmpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
        addpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));
        
        hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',~envelope_do));
        addpath(genpath(OSLDIR));
    end
    
    hmm.MixingMatrix = MixingMatrix;
    hmm.fsample = fsample;

    
    % Save results
    disp(['Saving inferred HMM to ' filenames.hmm])    
    save(filenames.hmm,'hmm');
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       O U T P U T   R E S U L T S                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.output
            
    statemaps = [filenames.output,'_',output_method];
 
    load(filenames.hmm);
    load(filenames.concat)

    switch output_method
        case {'pcorr'}
            
            D  = spm_eeg_load(filenames.prepare{1}); % to get number of voxels
            stat = zeros(D.nchannels,hmm.K);

            if use_parcels
                Dp = spm_eeg_load(prefix(filenames.prepare{1},'p'));
                statp = zeros(Dp.nchannels,hmm.K);

            end
                        
            epoched_statepath_sub = cell(length(data_files),1);
                        
            for subnum = 1:length(data_files)

                disp(['Computing ' output_method ' maps for ' data_files{subnum}]);
                
                % compute subject's state maps
                hmm_sub = hmm; 
                hmm_sub.statepath = hmm.statepath(subj_inds==subnum);                            
                hmm_sub = rmfield(hmm_sub,'MixingMatrix');
                
                
                D = spm_eeg_load(filenames.prepare{subnum});
                data = prepare_data(D,normalisation,logtrans);
                stat  = stat + osl_hmm_statemaps(hmm_sub,data,~envelope_do,output_method);
                
                if use_parcels
                    Dp = spm_eeg_load(prefix(filenames.prepare{subnum},'p'));
                    datap = prepare_data(D,normalisation,logtrans);
                    statp  = statp + osl_hmm_statemaps(hmm_sub,datap,~envelope_do,output_method);
                end
                
                good_samples = ~all(badsamples(D,':',':',':'));
                good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);

                sp_full = zeros(1,D.nsamples*D.ntrials);
                sp_full(good_samples) = hmm_sub.statepath;
                epoched_statepath_sub{subnum} = reshape(sp_full,[1,D.nsamples,D.ntrials]);
 
            end                        

            stat = stat ./ length(data_files);
            disp(['Saving state spatial maps to ' statemaps]) 
            
            nii_quicksave(stat,statemaps,getmasksize(D.nchannels),2);

                            
            if use_parcels
                % convert parcel statemaps into voxel statemaps
                pathstr = fileparts([filenames.prepare{1}]);
                fname   = fullfile(pathstr,'parcellation');
                load(fname); % Loads parcelAssignments and parcelWeights
                
                masksize = getmasksize(size(parcelAssignments,1));
                
                if ~use_parcel_weights
                    ROInets.nii_parcel_quicksave(statp,parcelAssignments,[statemaps,'_parcels'],masksize,2,'nearestneighbour');
                else % Question: Isn't it cheating to use the parcel weights to generate the maps?
                    weights = abs(parcelWeights);
                    weights = weights/mean(weights(logical(weights)));
                    
                    ROInets.nii_parcel_quicksave(statp,weights,[statemaps,'_parcels'],masksize,2);
                end
            end
                        
            % resave updated hmm
            hmm.statemaps = statemaps;
            hmm.epoched_statepath_sub = epoched_statepath_sub;
            disp(['Saving updated hmm to ' filenames.hmm]);
            save(filenames.hmm,'hmm');        

        otherwise
            warning([output_method ' is not a supported output method']);
    end

end


HMMresults = filenames.hmm;

end


function [data,cov_data] = prepare_data(D,normalisation,logtrans)

% reshape trialwise data
data = D(:,:,:);
data = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);

% select only good data
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
data = data(:,good_samples);

% log transform
if logtrans
    data = log10(data);
end

% normalisation
switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
end

cov_data = data;

end

