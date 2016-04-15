function [HMMresults,statemaps] = osl_hmm_groupinference_parcels(data_files,hmmdir,todo,options)
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
%              .prepare  - preps data, including computes amplitude prepares from source
%                          reconstructed data
%              .concat   - concatenates data and performs PCA dimensionality
%                          reduction
%              .hmm      - infers the HMM
%              .output   - output the state maps
%
%
% options    - struct with settings for preparing data and hmm stages
%              with the following fields:
%              .prepare - settings for prepare computation with fields:
%                          .envelope   - compute Hilbert envelope for data
%                                        (default 1)
%                          .freqbands  - frequencies within which to compute 
%                                        envelope {[Hz Hz],[Hz Hz]} only 
%                                        used if data is enveloped 
%                                        (default {} - wideband)
%                          .windowsize - window size for moving average filter [s]
%                                        only carried out if data is enveloped
%                                        (default 0.1)
%                          .filenames  - filenames to save/load prepares from
%                                        (default <hmmdir>/prepare/<BFfiles{subnum}>.mat)
%                          .parcellation.file - 3D niftii file with same
%                                               gridstep as data_files containing parcellation
%                                               labels. If not provided then no parcellation is
%                                               carried out (default)
%                          .parcellation.method - passed to ROInets.get_node_tcs
%                          .parcellation.protocol - passed to ROInets.remove_source_leakage
%                          .log - apply log transform to envelopes
%                          .mask_fname - mask used to get data from vols2matrix
%                                        [default is to assume a whole brain std space MNI mask]
%              .concat   - settings for temporal concatenation with fields:                      
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
%                          .nstates    - number of states to infer 
%                                        (default 8)
%                          .nreps      - number of repeat inferences to run 
%                                        (default 5)
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
    data_files = cellstr(data_files);
catch
    data_files=cellfun(@fullfile,data_files,'UniformOutput',false);
end;

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
try freqbands   = options.prepare.freqbands;  catch, freqbands   = {};  end
try logtrans    = options.prepare.log;        catch, logtrans    = 0;   end
try windowsize  = options.prepare.windowsize; catch, windowsize  = 0.1; end
try 
    options.prepare = ft_checkopt(options.prepare,'mask_fname','char');
catch
    options.prepare.mask_fname=[];
end;
mask_fname=options.prepare.mask_fname;

parcellation = [];
try parcellation.file     = options.prepare.parcellation.file;     catch, parcellation.file     = [];          end
try parcellation.protocol = options.prepare.parcellation.protocol; catch, parcellation.protocol = 'symmetric'; end
try parcellation.method   = options.prepare.parcellation.method;   catch, parcellation.method   = 'PCA';       end
try parcellation.normalise_voxeldata   = options.prepare.parcellation.normalise_voxeldata;   catch, parcellation.normalise_voxeldata   = 0;       end

use_parcels               = ~isempty(parcellation.file);

% Default concatenation settings
try pcadim        = options.concat.pcadim;        catch, pcadim         = 40;       end
try whiten        = options.concat.whiten;        catch, whiten         = 1;        end
try normalisation = options.concat.normalisation; catch, normalisation  = 'global'; end
try savePCmaps    = options.concat.savePCmaps;    catch, savePCmaps     = 0;        end

embed=[];
try embed.do      = options.concat.embed.do; catch, embed.do  = 0; end
try embed.centre_freq = options.concat.embed.centre_freq; catch, embed.centre_freq  = 15; end %Hz
try embed.rectify = options.concat.embed.rectify; catch, embed.rectify  = false; end %Hz

% Default HMM settings
try nstates     = options.hmm.nstates;      catch, nstates     = 8; end
try nreps       = options.hmm.nreps;        catch, nreps       = 5; end
try use_old_tbx = options.hmm.use_old_hmm_tbx;  catch, use_old_tbx = 0; end
try hmm_voxelwise = options.hmm.voxelwise;  catch, hmm_voxelwise = false; end

% Default output settings
try output_method       = options.output.method;                catch, output_method        = 'pcorr'; end
try state_assignment    = options.output.assignment;            catch, state_assignment     = 'hard';  end  
try use_parcel_weights  = options.output.use_parcel_weights;    catch, use_parcel_weights   = 0;       end


hmm=[];
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P   D A T A
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
if todo.prepare

    % check first subject to see if already parcellated
    D=spm_eeg_load(data_files{1});
    is_parcellated=isfield(D,'parcellation');
    if is_parcellated
        use_parcels=true;
    end;

    if is_parcellated
        parcelWeights=D.parcellation.weights;
        parcelAssignments=D.parcellation.assignments;
        try mask_fname=D.parcellation.mask_fname; catch mask_fname=[]; end;

        % Save parcellation results:
        pathstr = fileparts([filenames.prepare{1}]);
        fname = fullfile(pathstr,'parcellation');
        save(fname,'parcelAssignments','parcelWeights');
    end;
        
    for subnum = 1:length(data_files)                
        
        % Compute envelopes for all voxels
        if envelope_do
            disp('Computing Hilbert Envelopes')
        
            S         = [];
            S.D       = data_files{subnum};
            S.winsize = windowsize;
            S.freqbands = freqbands;
            D = osl_hilbenv(S);
        else
            % Need to copy D to preserve the original
            S=[];
            S.D=data_files{subnum};
            [pathstr,filestr] = fileparts(filenames.concat);
            S.outfile=[pathstr '/raw_' filestr];
            D = spm_eeg_copy(S);
        end
                
        % Move MEEG object to HMM directory
        move(D,filenames.prepare{subnum}); 
        disp(['Saving envelope data for ' D.fname ' to ' filenames.prepare{subnum}])
        clear D
             
        % Apply parcellation
        if use_parcels && ~is_parcellated

            disp('Applying parcellation')

            S                   = [];
            S.D                 = data_files{subnum};
            S.parcellation      = parcellation.file;
            S.orthogonalisation = parcellation.protocol;
            S.method            = parcellation.method;
            S.prefix            = 'p'; 
            S.normalise_voxeldata = parcellation.normalise_voxeldata;
            [D,parcelWeights,parcelAssignments] = osl_apply_parcellation(S);

            if subnum==1
                % Save parcellation results:
                pathstr = fileparts([filenames.prepare{1}]);
                fname = fullfile(pathstr,'parcellation');
                save(fname,'parcelAssignments','parcelWeights');
            end;
                    
            % Compute envelopes for all parcels
            if envelope_do
                S           = [];
                S.D         = fullfile(D.path,D.fname);
                S.winsize   = windowsize;
                S.freqbands = freqbands;
                D = osl_spmfun(@osl_hilbenv,S); % keep same filename
            end

            D = move(D,[fileparts(filenames.prepare{1}),filesep]);

        end

    end      

end

% for subnum = 1:length(data_files)
%     D=spm_eeg_load(filenames.prepare{subnum});
%     D.parcellation.S=S;
%     D.save;
% end;

% load parcel info
if use_parcels
    pathstr = fileparts([filenames.prepare{1}]);
    fname   = fullfile(pathstr,'parcellation');
    load(fname); % Loads parcelAssignments and parcelWeights
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    C O N C A T   H M M   D A T A                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.concat || (todo.infer && ~exist(filenames.concat,'file'))
            
    MixingMatrix = [];
    hmmdata      = [];
    
    for f = 1:max([numel(freqbands),1])
        % Load subjects and concatenate:
        dat_concat  = [];
        subj_inds   = [];
        C           = 0;
        Csamples    = 0;
        
        for subnum = 1:length(data_files)
            
            if use_parcels
                try 
                    D = spm_eeg_load(prefix(filenames.prepare{subnum},'p'));
                catch
                    warning('Loading in SPM MEG object that might not be parcellated.');
                    D = spm_eeg_load(filenames.prepare{subnum});
                end;
            else
                D = spm_eeg_load(filenames.prepare{subnum});
            end
            
            embed.tres=1/D.fsample;
                            
            %data = prepare_data(D,normalisation,logtrans,embed,'none');
            data = prepare_data(D,normalisation,logtrans,f,embed);
            
            % compute covariance (assuming zero mean data then global (over
            % all sessions) covariance matrix is:
            % (N1*C1 + N2*C2 + ... + Nk*Ck) / (N1 + N2 + ... + Nk)
            C = C + ((size(data,2)-1) * osl_cov(data));
            
            % keep track of number of samples
            Csamples = Csamples + size(data,2);
            
            dat_concat = [dat_concat, data];
            subj_inds  = [subj_inds, subnum*ones(1,size(data,2))];
            
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
        
        % Concatenate multiple frequency bands and mixing matrices
        hmmdata = [hmmdata (M * dat_concat)'];
        MixingMatrix = blkdiag(MixingMatrix,M);
        
        
        
        % Save PCA maps
        if ~embed.do && savePCmaps
            [pathstr,filestr] = fileparts(filenames.concat);
            if ~isempty(freqbands)
                bandstr = ['_fband' num2str(f)];
            else
                bandstr = '';
            end
            
            savefile_pc_maps        = [pathstr '/' filestr '_pc_maps'          bandstr];
            savefile_mean_pc_maps   = [pathstr '/' filestr '_mean_pc_maps'     bandstr];
            savefile_std_maps       = [pathstr '/' filestr '_hmmdata_std_maps' bandstr];
        
            if use_parcels                
                nii_settings            = [];
                nii_settings.interp     = 'nearestneighbour';
                nii_settings.mask_fname = mask_fname;
                ROInets.nii_parcel_quicksave(M', parcelAssignments, savefile_pc_maps,      nii_settings);
                ROInets.nii_parcel_quicksave(mean(abs(M),1)', parcelAssignments, savefile_mean_pc_maps, nii_settings);
                ROInets.nii_parcel_quicksave(std(pinv(M)*M*dat_concat,[],2), parcelAssignments, savefile_std_maps,     nii_settings);
            else      
                nii_settings            = [];
                nii_settings.mask_fname = mask_fname;
                nii_quicksave(M', savefile_pc_maps, nii_settings);
                nii_quicksave(mean(abs(M),1)', savefile_mean_pc_maps, nii_settings);
                nii_quicksave(std(pinv(M)*M*dat_concat,[],2), savefile_std_maps, nii_settings);
            end
            
            disp(['PCA maps saved to '          savefile_pc_maps]);
            disp(['mean of PCA maps saved to '  savefile_mean_pc_maps]);
            disp(['HMM data sd maps saved to '  savefile_std_maps]);
        end
          
        
    end %f
    
    fsample = D.fsample;
   
    disp(['Saving group concatenated data to ' filenames.concat])
    save(filenames.concat,'hmmdata','MixingMatrix','fsample','subj_inds','-v7.3');   
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            I N F E R   H M M                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.infer
    
    load(filenames.concat);

    if hmm_voxelwise
        numvox=size(D,1);
        numembeddings=size(hmmdata,2)/size(D,1);
        for vox=1:numvox
            start=(vox-1)*numembeddings+1;
            hmmdatavox=hmmdata(:,start:start+numembeddings-1);
            hmmdatavox=normalise(hmmdatavox,1);
            
            % Apply PCA dimensionality reduction            
            pcadim = floor(numembeddings/2);
            [allsvd,M] = eigdec(hmmdatavox'*hmmdatavox,pcadim);
       
            M = diag(1./sqrt(allsvd)) * M';
            hmmdatavox = [ (M * hmmdatavox')'];
        
            % Switch between Iead's & Diego's HMM toolboxes
            if use_old_tbx
                rmpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
                addpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

                if envelope_do || embed.rectify
                    hmm.hmm{vox} = ABhmm_infer(hmmdatavox,nstates,nreps);
                else
                    hmm.hmm{vox} = ABhmm_infer(hmmdatavox,nstates,nreps,'constrain_mean');
                end;
                addpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
                rmpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

            else

                rmpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
                addpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));

                hmm.hmm{vox} = osl_hmm_infer(hmmdatavox,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',~envelope_do));
                addpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
                rmpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));

            end
        end;
    else
        % Switch between Iead's & Diego's HMM toolboxes
        if use_old_tbx
            rmpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
            addpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

            if envelope_do || embed.rectify
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

    end;

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
        
    load(filenames.hmm);
    load(filenames.concat)
    
    switch output_method
        case {'pcorr','conn'}
            
            epoched_statepath_sub = cell(length(data_files),1);
             
            % hmm.statepath   = ABhmm_statepath(hmm);
            
            for f = 1:max([numel(freqbands),1])
                
                D  = spm_eeg_load(filenames.prepare{1}); % to get number of voxels
                stat = zeros(D.nchannels,hmm.K); 

                if ~isempty(freqbands)
                    bandstr = ['_fband' num2str(f)];
                else
                    bandstr = '';
                end
                
                statemaps{f} = [filenames.output,'_',output_method,bandstr];
                
                if use_parcels
                    try 
                        Dp = spm_eeg_load(prefix(filenames.prepare{1},'p'));
                    catch
                        warning('Loading in SPM MEG object that might not be parcellated.');
                        Dp = spm_eeg_load(filenames.prepare{1});
                    end;
                    statp = zeros(Dp.nchannels,hmm.K);
                end
                               
                % definitely do not want to do any embedding when computing pcorr
                % maps:
                embed.do=0;
                
                for subnum = 1:length(data_files)

                    disp(['Computing ' output_method ' maps for ' data_files{subnum}]);

                    % compute subject's state maps
                    hmm_sub = hmm;
                    hmm_sub.statepath = hmm.statepath(subj_inds==subnum);

%                    hmm_sub.train.Gamma=hmm.train.Gamma(subj_inds==subnum,:);
                    
                    hmm_sub = rmfield(hmm_sub,'MixingMatrix');

                    D = spm_eeg_load(filenames.prepare{subnum});

                    data = prepare_data(D,normalisation,logtrans,f,embed);
                    stat   = stat + osl_hmm_statemaps(hmm_sub,data,~envelope_do,output_method,state_assignment);
                    
                    if use_parcels
                        try 
                            Dp = spm_eeg_load(prefix(filenames.prepare{subnum},'p'));
                        catch
                            warning('Loading in SPM MEG object that might not be parcellated.');
                            Dp = spm_eeg_load(filenames.prepare{subnum});
                        end;
                        datap = prepare_data(Dp,normalisation,logtrans,f,embed);
                        statp  = statp + osl_hmm_statemaps(hmm_sub,datap,~envelope_do,output_method,state_assignment);
                    end

                    good_samples = ~all(badsamples(D,':',':',':'));
%                     good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);

                    if f==1
                        sp_full = zeros(1,D.nsamples*D.ntrials);
                        sp_full(good_samples) = hmm_sub.statepath;
                        epoched_statepath_sub{subnum} = reshape(sp_full,[1,D.nsamples,D.ntrials]);
                    end;
                end

                stat  = stat  ./ length(data_files);

                try
                    S2=[];
                    S2.mask_fname=mask_fname;
                    S2.output_spat_res=2; %mm
                    if isfield(D.parcellation.S,'hcp_sourcemodel3d') 
                        hcp_nii_quicksave(stat, [statemaps{f},'_wholebrain'], D.parcellation.S.hcp_sourcemodel3d, S2.output_spat_res, S2)
                    else
                        nii_quicksave(stat,[statemaps{f},'_wholebrain'],S2);
                    end;
                catch
                end;
                
                if use_parcels
                    statp = statp ./ length(data_files);
               

                    % convert parcel statemaps into voxel statemaps                    
                    S2.interp='nearestneighbour';

                    if isfield(D.parcellation,'S') && isfield(D.parcellation.S,'hcp_sourcemodel3d') 
                        hcp_nii_parcel_quicksave(statp, parcelAssignments, [statemaps{f},'_parcels'], D.parcellation.S.hcp_sourcemodel3d, S2.output_spat_res, S2)
                    else
                        ROInets.nii_parcel_quicksave(statp,parcelAssignments,[statemaps{f},'_parcels'],S2);
                    end
                    
                    % also store statemaps as vectors
                    hmm.statemap_parcel_vectors=statp;
                end
                
            end
            
            
            hmm.epoched_statepath_sub = epoched_statepath_sub;
            hmm.statemaps = statemaps;
            hmm.filenames=filenames;
                         
        otherwise
            
            warning([output_method ' is not a supported output method']);
    end
    
end

% save updated hmm
hmm.filenames=filenames;
disp(['Saving updated hmm to ' filenames.hmm]);
save(filenames.hmm,'hmm');
           
HMMresults = filenames.hmm;

end

function data = prepare_data(D,normalisation,logtrans,freq_ind,embed)

% reshape trialwise data
if strcmp(D.transformtype,'TF')
    data = D(:,freq_ind,:,:);
else
    data = D(:,:,:);
end
data = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);

% select only good data
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
data = data(:,good_samples);

% log transform
if logtrans
    data = log10(data);
end

if exist('embed','var') && embed.do,    
    
    disp('Time embedding data');
    span=1/embed.centre_freq; %secs
    num_embeddings=round(span/embed.tres);        
    %lags=round(linspace(-num_embeddings/2,num_embeddings/2,num_embeddings));
    lags=round(-num_embeddings/2:num_embeddings/2);
    %lags=round(0:num_embeddings-1);

    disp(lags);
    
    [dataout,valid]=embedx(data',lags);
    dataout=dataout';
    data=randn(size(dataout,1),size(data,2))*std(squash(data(:,1:500)));
    data(:,valid)=dataout;          
       
end;

% normalisation
switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
end

if embed.rectify
    data=abs(data);
    
    switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
    end;
end

end

