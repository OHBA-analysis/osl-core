function [HMMresults,statemaps] = teh_groupinference_parcels(data_files,options)
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
%              .prepare  - preps data, including performs PCA dimensionality
%                          reduction
%              .hmm      - infers the HMM
%              .output   - output the state maps
%
%
% options    - struct with settings for preparing data and hmm stages
%              with the following fields:
%              
%              .prepare   - settings for prep with fields:                      
%                          .pcadim        - dimensionality to use (default 40)             
%                          .normalisation - variance normalisation method:
%                            'none'     [data = demean(data,2)]
%                            'global'   [data = demean(data,2)/std(data(:)]
%                            'voxelwise'[data = normalise(data,2)]
%                            (default 'global')
%                          .whiten        - apply whitening [0/1] (default 1)
%                          .embed.do      - do time embedding using call to
%                                           embedx 
%                          .embed.centre_freq - centre freq of interest (Hz), used to set the time span of the time embedding 
%                                               (default is 15Hz)
%                                               (default is 0)
%                          .filename      - filename to save/load concatenated data from
%                                           (default <hmmdir>/env_concat.mat)
%                          .savePCmaps    - save PCA maps [0/1] (default 0)
%
%              .hmm      - settings for HMM inference with fields:                
%                          .options    - options to pass to hmmmar
%                          .filename   - filename to save/load HMM results from
%                                        (default <hmmdir>/HMM.mat)
%      
%              .output   - settings for the HMM output with fields
%                          .method     - maps/matrix to output:
%                                        'pcorr' - partial correlation maps           
%                                        (default 'pcorr')
%                          .assignment - Use hard or soft (probabilistic)
%                                        state assignment ['soft','hard']
%                                        (default 'hard')
%                          .use_parcel_weights - uses parcel weights rather
%                          than binary parcels [0/1] (default 1)
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
% Mostly AB & a little bit from MWW (2014/15)

global OSLDIR

statemaps=[];

try 
    todo = options.todo;
catch
    todo=[];
    todo.prepare=1;
    todo.hmm=1;
    todo.output=1;
end

try 
    hmmdir = options.hmmdir;
catch 
    error('options.hmmdir must be specified');
end

try 
    data_files = cellstr(data_files);
catch
    data_files=cellfun(@fullfile,data_files,'UniformOutput',false);
end

data_fnames = cell(size(data_files));
for f = 1:numel(data_files)
  [~,data_fnames{f},~] = fileparts(data_files{f});
end

if ~isdir(hmmdir)
  mkdir(hmmdir); 
end

if isfield(options.prepare,'filename') && ~isempty(options.prepare.filename)
    [pathstr,filestr] = fileparts(options.prepare.filename);
        if strcmp(pathstr,'')
            pathstr = hmmdir;
        end
        filenames.prepare = fullfile(pathstr,filestr);
else
    filenames.prepare = fullfile(hmmdir,'prepare');
end
filenames.prepare = [filenames.prepare '.mat'];

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
try todo.hmm     = todo.hmm;     catch, todo.hmm     = 1; end

% Default prepare settings
try pcadim        = options.prepare.pcadim;        catch, pcadim         = 40;       end
try whiten        = options.prepare.whiten;        catch, whiten         = 1;        end
try normalisation = options.prepare.normalisation; catch, normalisation  = 'global'; end
try savePCmaps    = options.prepare.savePCmaps;    catch, savePCmaps     = 0;        end
logtrans=0;

embed=[];
try embed.do      = options.prepare.embed.do; catch, embed.do  = 0; end
try embed.centre_freq = options.prepare.embed.centre_freq; catch, embed.centre_freq  = 15; end %Hz
try embed.rectify = options.prepare.embed.rectify; catch, embed.rectify  = false; end %Hz

% Default HMM settings
try hmmoptions = options.hmm;  catch, hmmoptions=[]; end

% Default output settings
try output_method       = options.output.method;                catch, output_method        = 'pcorr'; end
try state_assignment    = options.output.assignment;            catch, state_assignment     = 'hard';  end  
try use_parcel_weights  = options.output.use_parcel_weights;    catch, use_parcel_weights   = 0;       end

hmm=[];
hmm.options=options;

%%%%%%%%%%%%%%%%%%%%%
%% check parcellation and extract parcellation info

% check first subject to see if it is parcellated
D=spm_eeg_load(data_files{1});
is_parcellated=isfield(D,'parcellation');
if ~is_parcellated
    error('Data must be parcellated');
end;

parcelWeights=D.parcellation.weights;
parcelAssignments=D.parcellation.assignments;
try mask_fname=D.parcellation.mask_fname; catch mask_fname=[]; end

%%%%%%%%%%%%%%%%%%%%%
%% create filenames for prepared data for each subj
%filenames.prepared_data=[];
%[pathstr] = fileparts(filenames.prepare);
%for subnum = 1:length(data_files)
%    [~,filestr] = fileparts(data_files{subnum});
%    filenames.prepared_data{subnum}=[pathstr '/prepared_' filestr '.mat'];
%end

[pathstr] = fileparts(filenames.prepare);
filenames.prepared_data=[pathstr '/prepared_data.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P A R E   H M M   D A T A                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.prepare
            
    %MixingMatrix = [];
    hmmdata      = [];

    % Load subjects and concatenate:
    dat_concat  = [];
    subj_inds   = [];
    C           = 0;
    Csamples    = 0;
    hmmT=[];
    
    for subnum = 1:length(data_files)

        D = spm_eeg_load(data_files{subnum});            

        embed.tres=1/D.fsample;

        data = osl_teh_prepare_data(D,normalisation,logtrans,f,embed);

        % compute covariance (assuming zero mean data then global (over
        % all sessions) covariance matrix is:
        % (N1*C1 + N2*C2 + ... + Nk*Ck) / (N1 + N2 + ... + Nk)
        C = C + ((size(data,2)-1) * osl_cov(data));

        % keep track of number of samples
        Csamples = Csamples + size(data,2);

        dat_concat = [dat_concat, data];
        subj_inds  = [subj_inds, subnum*ones(1,size(data,2))];
        
        hmmT{subnum} = length(data);

    end

    C = C / (Csamples-1);

    clear data;

    % Apply PCA dimensionality reduction
    if pcadim > 0
        disp(['Keeping top ' num2str(pcadim) ' PCs']);
        pcadim = min(pcadim,size(dat_concat,1));
        [allsvd,M] = eigdec(C,pcadim);
    else
        M = eye(size(dat_concat,1));
        allsvd = diag(C);
    end

    if whiten
        M = diag(1./sqrt(allsvd)) * M';
    else
        M = M';
    end

    % Apply PC dim to concat data
    hmmdata = [hmmdata (M * dat_concat)'];
    %MixingMatrix = blkdiag(MixingMatrix,M);
    
    % Apply PC dim to each subj separately 
    clear X;
    for subnum = 1:length(data_files)

        D = spm_eeg_load(data_files{subnum});            

        embed.tres=1/D.fsample;

        [data, lags] = osl_teh_prepare_data(D,normalisation,logtrans,f,embed);
    
        X{subnum}=(M * data)';          
    end
    
    if 0
        
        %% 
        %remove middle block:
        C2=triu(C,length(lags))+tril(C,-length(lags));
        MC2M=pinv(M)*M*C*(pinv(M)*M)';
        MC2M2=triu(MC2M,length(lags))+tril(MC2M,-length(lags));
        figure;subplot(121);imagesc(C2, [-.1 .1]);colorbar;subplot(122);imagesc(MC2M2, [-.1 .1]);colorbar;
        figure;subplot(121);hist(squash(C2(C2~=0)),20);subplot(122);hist(squash(MC2M2(MC2M2~=0)),20);
        % middle block:
        C3=C-C2;
        MC2M3=MC2M-MC2M2;
        figure;subplot(121);hist(squash(C3(C3~=0)));subplot(122);hist(squash(MC2M3(MC2M3~=0)),20);
    end
    
    % Save PCA maps
    if savePCmaps        
        if ~embed.do
            [pathstr,filestr] = fileparts(filenames.prepare);

            savefile_pc_maps        = [pathstr '/' filestr '_pc_maps'          ];
            savefile_mean_pc_maps   = [pathstr '/' filestr '_mean_pc_maps'     ];
            savefile_std_maps       = [pathstr '/' filestr '_hmmdata_std_maps' ];

            nii_settings            = [];
            nii_settings.interp     = 'nearestneighbour';
            nii_settings.mask_fname = mask_fname;
            ROInets.nii_parcel_quicksave(M', parcelAssignments, savefile_pc_maps,      nii_settings);
            ROInets.nii_parcel_quicksave(mean(abs(M),1)', parcelAssignments, savefile_mean_pc_maps, nii_settings);
            ROInets.nii_parcel_quicksave(std(pinv(M)*M*dat_concat,[],2), parcelAssignments, savefile_std_maps,     nii_settings);

            disp(['PCA maps saved to '          savefile_pc_maps]);
            disp(['mean of PCA maps saved to '  savefile_mean_pc_maps]);
            disp(['HMM data sd maps saved to '  savefile_std_maps]);
        end
    end
   
    fsample = D.fsample;
   
    disp(['Saving prepared data to ' filenames.prepare])
    %save(filenames.concat,'hmmdata','MixingMatrix','fsample','subj_inds','-v7.3');   
    save(filenames.prepare,'lags','M','C','hmmT','hmmdata','fsample','subj_inds','-v7.3');   
    save(filenames.prepared_data,'X','-v7.3');
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            I N F E R   H M M                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.hmm
           
    % load in what's needed
    if ~todo.prepare
        load(filenames.prepare,'hmmT','hmmdata','fsample','subj_inds');  
        load(filenames.prepared_data,'X');  
    end
    
    % setup default options for hmm
    options = [];
    options.inittype = 'hmmmar';
    options.initcyc = 30;
    options.initrep = 1;
    options.order = 0;
    options.zeromean = 1;
    options.K=8;

    if hmmoptions.big
        % defaults for stochastic hmm
        options.BIGNbatch = 5;
        options.BIGtol = 1e-7;
        options.BIGcyc = 500;
        options.BIGundertol_tostop = 5;
        options.BIGdelay = 5; 
        options.BIGforgetrate = 0.7;
        options.BIGbase_weights = 0.9;
        options.BIGuniqueTrans = 1;
    end
    
    % copy passed in hmmoptions over default options 
    fn = fieldnames(hmmoptions);
    for fi=1:numel(fn)
        options.(fn{fi}) = hmmoptions.(fn{fi});
    end
    % remove options not used by hmmmar call
    options=rmfield(options,'big');
    if ~hmmoptions.big && isfield(options,'BIGNbatch')
        options=rmfield(options,'BIGNbatch');        
    end    
    
    % call hmmmar
    [hmm, ~, ~, ~, ~, ~, fehist] = hmmmar(X,hmmT,options);
    
    % save('/Users/woolrich/for_diego','X','hmmT','options'); !scp /Users/woolrich/for_diego.mat $WS12HD
    % Error running [hmm, ~, ~, ~, ~, ~, fehist] = hmmmar(X,hmmT,options); You can get the data from: hbaws12.ohba.ox.ac.uk:/home/mwoolrich/homedir/for_diego.mat
    
    hmm.options = hmmoptions;
    hmm.data_files = data_files;
    hmm.fsample = fsample;
    hmm.fehist = fehist;
    %hmm.gamma = gamma;
    %hmm.statepath = vpath;
    
    disp(['Saving inferred HMM to ' filenames.hmm])    
    save(filenames.hmm,'hmm');

    % compute state time courses    
    hmm.gamma = hmmdecode(X,hmmT,hmm,0);  % last argument: 0, state time courses; 1, viterbi path
    hmm.statepath = hmmdecode(X,hmmT,hmm,1);  % last argument: 0, state time courses; 1, viterbi path
    
    %hmm.statepath=zeros(size(hmm.statepath_hot,1),1);
    %for ii=1:size(hmm.statepath_hot,1)
    %    ind = find(hmm.statepath_hot(ii,:));
    %    hmm.statepath(ii)=ind;
    %end;
    
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
    
    % load in what's needed
    if ~todo.prepare
        load(filenames.prepare);  
    end
    if ~todo.hmm
        load(filenames.hmm);
    end

    switch output_method
        case {'pcorr','conn'}
            
            epoched_statepath_sub = cell(length(data_files),1);
             
            statemaps = [filenames.output,'_',output_method];

            % load first subj to check 
            Dp = spm_eeg_load(data_files{1});
            statp = zeros(Dp.nchannels,hmm.K);

            % definitely do not want to do any embedding when computing pcorr
            % maps:
            embed.do=0;

            hmm.statemap_parcel_vectors_persubj=zeros(length(data_files),Dp.nchannels,hmm.K);
            for subnum = 1:length(data_files)

                disp(['Computing ' output_method ' maps for ' data_files{subnum}]);

                % compute subject's state maps
                hmm_sub = hmm;
                %hmm_sub.gamma = hmm.gamma(subj_inds==subnum,:);
                hmm_sub.statepath = hmm.statepath(subj_inds==subnum);
                                           
                Dp = spm_eeg_load(data_files{subnum});
                datap = osl_teh_prepare_data(Dp,normalisation,logtrans,f,embed);
                
                tmp=osl_hmm_statemaps(hmm_sub,datap,true,output_method,state_assignment);            

                statp  = statp + tmp;
                hmm.statemap_parcel_vectors_persubj(subnum,:,:)=tmp;
                
                good_samples = ~all(badsamples(Dp,':',':',':'));

                sp_full = zeros(1,Dp.nsamples*Dp.ntrials);
                sp_full(good_samples) = hmm_sub.statepath;
                epoched_statepath_sub{subnum} = reshape(sp_full,[1,Dp.nsamples,Dp.ntrials]);
                
                if Dp.ntrials>1
                    hmm.is_epoched=1;
                else
                    hmm.is_epoched=0;
                end
            end

            statp = statp ./ length(data_files);

            % convert parcel statemaps into voxel statemaps                    
            S2.interp='nearestneighbour';

            if isfield(D.parcellation,'S') && isfield(D.parcellation.S,'hcp_sourcemodel3d') 
                hcp_nii_parcel_quicksave(statp, parcelAssignments, [statemaps,'_parcels'], D.parcellation.S.hcp_sourcemodel3d, S2.output_spat_res, S2)
            else
                ROInets.nii_parcel_quicksave(statp,parcelAssignments,[statemaps,'_parcels'],S2);
            end

            % also store statemaps as vectors
            hmm.statemap_parcel_vectors=statp;        
            
            hmm.epoched_statepath_sub = epoched_statepath_sub;
            hmm.statemaps = statemaps;
            hmm.filenames=filenames;

        otherwise
            
            warning([output_method ' is not a supported output method']);
            
    end
    
end

% save updated hmm
hmm.filenames=filenames;
hmm.data_files=data_files;

disp(['Saving updated hmm to ' filenames.hmm]);
save(filenames.hmm,'hmm');
           
HMMresults = filenames.hmm;

end


