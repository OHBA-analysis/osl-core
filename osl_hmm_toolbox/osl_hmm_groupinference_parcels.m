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
%                            (default <hmmdir>/prepares/<BFfiles{subnum}>.mat)
%                          .parcellation.file - 3D niftii file with same
%                          gridstep as data_files containing parcellation
%                          labels. If not provided then no parcellation is
%                          carried out (default)
%                          .parcellation.method - passed to ROInets.get_corrected_node_tcs
%                          .parcellation.protocol - passed to ROInets.get_corrected_node_tcs
%                          .log        - apply log transform to envelopes
%              .concat   - settings for temporal concatenation with fields:%                          
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
try todo.concat   = todo.concat;   catch, todo.concat   = 1; end
try todo.hmm      = todo.hmm;      catch, todo.hmm      = 1; end

% Default prepare settings
try windowsize = options.prepare.windowsize; catch, windowsize = 0.1;  end
parcellation=[];
try parcellation.file  = options.prepare.parcellation.file;  catch, parcellation.file  = [];   end
use_parcels = ~isempty(parcellation.file);
try parcellation.protocol  = options.prepare.parcellation.protocol;  catch, parcellation.protocol  = 'symmetric';   end
try parcellation.method  = options.prepare.parcellation.method;  catch, parcellation.method  = 'PCA';   end
try envelope_do  = options.prepare.envelope;  catch, envelope_do  = 1;   end
try logtrans  = options.prepare.log;        catch, logtrans   = 0;  end

% Default concatenation settings
try pcadim    = options.concat.pcadim;     catch, pcadim     = 40; end
try whiten    = options.concat.whiten;     catch, whiten     = 1;  end

% Default HMM settings
try nstates = options.hmm.nstates; catch, nstates = 8; end
try nreps   = options.hmm.nreps;   catch, nreps   = 5; end
    
% Default output settings
try output_method = options.output.method; catch, output_method = 'pcorr'; end
try use_parcel_weights = options.output.use_parcel_weights; catch, use_parcel_weights = 1; end
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P   D A T A
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.prepare

    if use_parcels, 
        
        % creat parcellation.parcelflag 
        parcels=read_avw(parcellation.file);
        D = spm_eeg_load(data_files{1});
        D = montage(D,'switch',2); D.save;  
        parcellation.gridstep=getmasksize(D.nchannels);
        stdbrain = read_avw([OSLDIR '/std_masks/MNI152_T1_' num2str(getmasksize(D.nchannels)) 'mm_brain.nii.gz']);  
        
        try
            parcellation.parcelflag=vols2matrix(parcels,stdbrain); %nVoxels x nSamples
        catch
            error('Make sure options.prepare.parcellation.file and the data are compatible, including having the same spatial resolution.');
        end;
        
        if size(parcellation.parcelflag,2)==1
            % Currently nvoxels x 1 with a index at each voxel indicating 
            % parcel membership.
            % Instead, needs to be nvoxels x nparcels binary labelling
            num_parcels=max(parcellation.parcelflag);
            parcelflag_new=zeros(size(parcellation.parcelflag,1),num_parcels);
            for indind=1:size(parcellation.parcelflag,1),
                if parcellation.parcelflag(indind)>0
                    parcelflag_new(indind,parcellation.parcelflag(indind))=1;
                end;
            end;
            parcellation.parcelflag=parcelflag_new;
            clear parcelflag_new;
        end;
        
        % create parcellation based on 
        % winner takes all voting (this will be used for outputting
        % results in nii format)
        parcellation.parcelflag_voted=logical(zeros(size(parcellation.parcelflag)));
        for indind=1:size(parcellation.parcelflag,1),
            [a b]=max(parcellation.parcelflag(indind,:));
            if a>0
               parcellation.parcelflag_voted(indind, b)=1;
            end;
        end;  
    else
        parcellation=[];
    end;    
  
    for subnum = 1:length(data_files)
    
        [pathstr,filestr] = fileparts(data_files{subnum});

        D = spm_eeg_load(data_files{subnum});
        D = montage(D,'switch',2); 
        D.save;

        %%% get voxel data   
        disp(['Computing voxel data for ' filestr]);
        S2=[];
        S2.D=D;
        S2.montage_index=2; % with weights normalisation        
        S2.D_block.size=min(1000,size(D,1));
        data=zeros(size(D,1),size(D,2),size(D,3));

        ft_progress('init', 'etf');
        for indind=1:size(D,1), % indexes brain space  
            ft_progress(indind/size(D,1));
            S2.index=indind;
            [data(indind,:,:) S2] = osl_get_recon_timecourse( S2 );                       
        end; 
        ft_progress('close');

        % remove bad trials and samples
        disp(['Removing bad trials and samples for ' filestr]);
        data=clean_data(D,data);        

        ntrials=size(data,3); % could be different for each subject

        if subnum==1, % these should be the same for every subject:
            nsamples=size(data,2);
            nvoxels=size(data,1);
        end;

        if use_parcels,  
            disp(['Computing node data for ' filestr]);

            %data2=devar(data,2);
            data2=data;
            
            voxeldata_concat = reshape(data2,[nvoxels,nsamples*ntrials]);
            
            % sds=std(voxeldata_concat,[],2);nii_quicksave(sds,'~/tmpvox',8,8);
            % !fslview ~/tmpvox &
            
            [nodedata_concat parcellation.voxelWeightings] = ROInets.get_corrected_node_tcs(voxeldata_concat, parcellation.parcelflag, parcellation.protocol, parcellation.method);

            % sds=std(nodedata_concat,[],2);ROInets.nii_parcel_quicksave(sds,parcellation.parcelflag_voted,'~/tmpparc',8,8);
            % !fslview ~/tmpparc &
 
            if subnum==1, % this should be the same for every subject:
                nnodes=size(nodedata_concat,1);
            end;

            data=reshape(nodedata_concat,[nnodes,nsamples,ntrials]);
        end;
        
        clear voxeldata_concat nodedata_concat;

        if logtrans 
            data = log10(data);
        end 

        % save raw data 
        [pth nm ext]=fileparts(filenames.prepare{subnum});
        outfile=[pth '/' nm '_raw' ext];    
        disp(['Saving raw data for ' filestr ' to ' outfile])
        Dclean = clone(montage(D,'switch',0),outfile,[size(data,1),nsamples,ntrials]);
        Dclean(:,:,:) = data;
        Dclean.save;

        if envelope_do
            disp(['Computing envelope data for ' filestr]);

            ft_progress('init', 'etf');                
            for trl = 1:ntrials,                 
                ft_progress(trl/ntrials);
                [tmp t_env] = hilbenv(data(:,:,trl),D.time,windowsize*D.fsample);                   
                if trl==1,
                    env_data=zeros(size(data,1), size(tmp,2), ntrials);
                end;
                env_data(:,:,trl)=tmp;              
            end;
            ft_progress('close'); 

            if logtrans 
                env_data = log10(env_data);
            end 

            % save env data
            [pth nm ext]=fileparts(filenames.prepare{subnum});
            outfile=[pth '/' nm '_env' ext];    
            disp(['Saving envelope data for ' filestr ' to ' outfile])
            Dclean = clone(montage(D,'switch',0),outfile,[size(env_data,1),size(env_data,2),ntrials]);
            Dclean(:,:,:) = env_data;
            Dclean = timeonset(Dclean,t_env(1));
            fsampleNew = 1./(diff(t_env([1,end]))/length(t_env));
            Dclean = Dclean.fsample(fsampleNew);
            Dclean.save;

        end;   

    end;
      
    [pth fn ext]=fileparts([filenames.prepare{1}]);
    fname=[pth '/' fn '_parcellation'];
    disp(['Saving parcellation to ' fname])
    save(fname,'parcellation');
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    C O N C A T   H M M   D A T A                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.concat || (todo.infer && ~exist(filenames.concat,'file'))
            
    % Load subjects and concatenate:
    dat_concat = [];
    C = 0;
    Csamps = 0;
    
    for subnum = 1:length(data_files)

        S=[];
        S.envelope_do=envelope_do;
        S.filename=filenames.prepare{subnum};                
        D=load_subject_data(S);        
        data = D(:,:,:);
        
        % normalisation
        if envelope_do
            data = demean(data,2)./std(data(:));
            cov_data=data;
            %cov_data = normalise(data,2);
        else
            data = normalise(data,2);
            cov_data = data;
        end;

        ntrials=size(data,3); % could be different for each subject

        if subnum==1, % these should be the same for every subject:
            nsamples=size(data,2);
            nvoxels=size(data,1);
            fsample=D.fsample;
        end;

        % compute covariance
        for trl=1:ntrials,                            
            cov_data_tri=cov_data(:,:,trl);
            C = C + cov_data_tri*permute(cov_data_tri,[2,1]); 
        end;        
        
        % concat data over trials
        data=reshape(data,[nvoxels,nsamples*ntrials]);
        Csamps=Csamps+nsamples*ntrials;
        
        dat_concat=[dat_concat, data];
    end
    
    C=C/Csamps;
    
    clear data cov_data;

    if pcadim>0,
        
        pcadim=min(pcadim,size(dat_concat,1));
            
        disp(['Keeping top ' num2str(pcadim) ' PCs']);
        
        % PCA + whitening - below is equivalent to fastica whitening code but much faster
        [allsvd,MixingMatrix] = eigdec(C,pcadim);
        
    else
        MixingMatrix = eye(size(dat_concat,1));   
        allsvd = diag(C);
    end;
    
    if whiten
        MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
    else
        MixingMatrix = MixingMatrix';
    end
    
    hmmdata = (MixingMatrix * dat_concat)';
            
    disp(['Saving group concatenated data to ' filenames.concat])
    save(filenames.concat,'hmmdata','MixingMatrix','fsample');

    %%%%%%%%%%
    % Save mixing matrix 

    disp(['Saving mixing matrix'])

    % get parcelflag (if needed) and masksize
    if use_parcels
        [pth fn ext]=fileparts([filenames.prepare{1}]);
        fname=[pth '/' fn '_parcellation'];
        load(fname,'parcellation');
        msize=parcellation.gridstep;
    else
        msize=getmasksize(size(dat_concat,1));
    end;

    [pth nm ext]=fileparts(filenames.concat);
    pc_maps=[pth '/' nm '_pc_maps'];
    
    if use_parcels
        pc_maps = ROInets.nii_parcel_quicksave(MixingMatrix',parcellation.parcelflag_voted,pc_maps,msize,msize);
    else
        pc_maps = nii_quicksave(MixingMatrix',pc_maps,msize,2);
    end;
    disp(['PCA maps saved to ' pc_maps]);   

    [pth nm ext]=fileparts(filenames.concat);
    mean_pc_maps=[pth '/' nm '_mean_pc_maps'];

    sd_maps=[pth '/' nm '_hmmdata_sd_maps'];
    sds = std(pinv(MixingMatrix)*hmmdata',[],2);

    if use_parcels
        pc_maps = ROInets.nii_parcel_quicksave(mean(abs(MixingMatrix),1)',parcellation.parcelflag_voted,mean_pc_maps,msize,2);
        sd_maps = ROInets.nii_parcel_quicksave(sds,parcellation.parcelflag_voted,sd_maps,msize,2);
    else
        pc_maps = nii_quicksave(mean(abs(MixingMatrix),1)',mean_pc_maps,msize,msize);
        sd_maps = nii_quicksave(sds,sd_maps,msize,msize);
    end;
    disp(['mean of PCA maps saved to ' mean_pc_maps]);  
    disp(['HMM data sd maps saved to ' sd_maps]);     

    %%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            I N F E R   H M M                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.infer
    
    load(filenames.concat);
    D = spm_eeg_load(data_files{1}); 
    
    %%%%%%%%%%%%%%%%
    % Sort this OUT!
    if options.hmm.use_old_hmm_tbx
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       O U T P U T   R E S U L T S                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.output
            
    statemaps = [filenames.output,'_',output_method];

    HMMresults = filenames.hmm;
    load(HMMresults)

    switch output_method
        case {'pcorr'}
            % Compute partial correlation/variance maps
            stat = zeros(size(hmm.MixingMatrix,2),hmm.K);                

            subjstart_index=1;                

            epoched_statepath_sub=cell(length(data_files),1);
                        
            for subnum = 1:length(data_files)

                disp(['Computing ' output_method ' maps for ' data_files{1}]);
                
                dat_concat_sub=[];                                                
            
                S=[];
                S.envelope_do=envelope_do;
                S.filename=filenames.prepare{subnum};
                D=load_subject_data(S);
                data=D(:,:,:);
                
                % normalisation
                if envelope_do
                    if pcadim>0,
                        data = demean(data,2)./std(data(:));                        
                    else
                        data = demean(data,2)./std(data(:));
                    end;
                else
                    data = normalise(data,2);
                end;

                ntrials=size(data,3); % could be different for each subject

                if subnum==1, % these should be the same for every subject:
                    nsamples=size(data,2);
                    nvoxels=size(data,1);
                end;

                % get concat data over trials
                dat_concat_sub=reshape(data,[nvoxels,nsamples*ntrials]);
                      
                % extract subject's state paths
                from=subjstart_index;
                to=from+size(dat_concat_sub,2)-1;
                subjstart_index=to+1;
                epoched_statepath_sub{subnum}=hmm.statepath(from:to);
            
                % compute subject's state maps
                hmm_sub = hmm; 
                hmm_sub.statepath = epoched_statepath_sub{subnum};                              
                hmm_sub=rmfield(hmm_sub,'MixingMatrix');
                
                stat = stat + osl_hmm_statemaps(hmm_sub,dat_concat_sub,~envelope_do,output_method);

                % unconcat epoched_statepath_sub
                epoched_statepath_sub{subnum}=reshape(epoched_statepath_sub{subnum},[1,nsamples,ntrials]);
            end                        

            stat = stat ./ length(data_files);
            disp(['Saving state spatial maps to ' statemaps])    
            if use_parcels
                % convert parcel statemaps into voxel statemaps
                [pth fn ext]=fileparts([filenames.prepare{1}]);
                fname=[pth '/' fn '_parcellation'];
                load(fname,'parcellation');
                msize=parcellation.gridstep;
                
                if ~use_parcel_weights
                    statemaps = ROInets.nii_parcel_quicksave(stat,parcellation.parcelflag_voted,statemaps,msize,2);
                else
                    weights=abs(parcellation.parcelflag);
                    weights=weights/mean(weights(logical(weights)));                    

                    statemaps = ROInets.nii_parcel_quicksave(stat,weights,statemaps,msize,2);
                end;
            else
                statemaps = nii_quicksave(stat,statemaps,getmasksize(D.nchannels),2);
            end;
                        
            % resave updated hmm
            hmm.statemaps=statemaps;
            hmm.epoched_statepath_sub=epoched_statepath_sub;
            disp(['Saving updated hmm to ' filenames.hmm]);
            save(filenames.hmm,'hmm');        

        otherwise
            warning([output_method ' is not a supported output method']);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clean_data=clean_data(D,data)
    
% extracts removes bad trials and segments

clean_data=[];

trials = D.indtrial(D.condlist,'good'); 

for trl=1:length(trials),                            

    tbad = all(badsamples(D,':',':',trials(trl)));
    samples2use = find(~tbad);
    clean_data(:,samples2use,trl) = data(:,samples2use,trials(trl)); %#ok - SPM doesn't like logical indexing        
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data=load_subject_data(S)

if S.envelope_do
    % load env data
    [pth nm ext]=fileparts(S.filename);
    outfile=[pth '/' nm '_env' ext];    
    disp(['Loading envelope data ' outfile])
else
    % load raw data
    [pth nm ext]=fileparts(S.filename);
    outfile=[pth '/' nm '_raw' ext];    
    disp(['Loading raw data ' outfile])         
end;

data = spm_eeg_load(outfile);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



