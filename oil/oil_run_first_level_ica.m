function  oil = oil_run_first_level_ica( oil )

% [ oil ] = osl_run_first_level_ica( oil )
%
% The function runs basic subject level continuous GLM analyses on the ICA
% outputs from osl_ica.
%
% This function should normally be called using osl_run_ica(oil);
%
% HL 05.02.13 version 1.1

OSLDIR = getenv('OSLDIR');
switch lower(oil.ica.temp_or_spat)
    case {'temporal'}
        if isfield(oil.ica.results,'tICs'); NumICs=size(oil.ica.results.tICs,1);end
    case {'spatial'}
        if isfield(oil.ica.results,'tICs'); NumICs=size(oil.ica.results.mixing_matrix,1);end
    otherwise
end

first_level=oil.ica_first_level;

contrast_list=first_level.contrast;
for c=1:length(contrast_list),
    contrast_list{c}=contrast_list{c}(:);
end;

clear stats;

if(~isfield(first_level,'design_matrix') && ~isfield(first_level,'design_matrix_summary')),
    error('Design matrix is not specified');
end;

for sub=1:length(oil.source_recon.results_fnames(oil.concat_subs.sessions_to_do)), % indexes subjects
        
    % load in the beamformer result for this subject
    source_recon_results_fname = oil.source_recon.results_fnames{oil.concat_subs.sessions_to_do(sub)};
    source_recon_results=osl_load_oil_results(oil,source_recon_results_fname);
    trialtypes=[];
    triallist=[];
    D=source_recon_results.BF.data.D;
    for i=1:length(source_recon_results.source_recon.conditions), % indexes conditions/triggers within subset
        Ntrialspercond=0;
        Ntrialspercond=Ntrialspercond+length(D.pickconditions(source_recon_results.source_recon.conditions{i})); %% number of trials for this condition
        triallist=[triallist , D.pickconditions(source_recon_results.source_recon.conditions{i})];
        trialtypes=[trialtypes ; ones(Ntrialspercond,1).*i];
    end

    Nsamples=D.nsamples;
    Ntrials=length(trialtypes);

    if(Ntrials==1), error('Number of trials should be more than 1 for a epoched analysis');   end;
    
    % Pre-allocate results containers
    if(sub==1),
        first_level_results.stdcope=zeros([NumICs,length(contrast_list)]);% num_components x num_contrasts
        first_level_results.cope=zeros([NumICs,length(contrast_list)]);% num_components x  num_contrasts
    end;
        
    % setup the GLM design matrix for this subject
    if(~isfield(first_level,'design_matrix')),
        S3=[];
        S3.trialtypes=trialtypes;
        S3.Xsummary=first_level.design_matrix_summary;
        if(isfield(first_level,'trial_rejects'))
            S3.trial_rejects=first_level.trial_rejects;
        end;
        first_level.x=setup_beamformer_designmatrix(S3);
    end;
    
    % precompute some things
    contrast_list=first_level.contrast;
    x=first_level.x;
    
    
    Nsamples_per_trial=floor((source_recon_results.source_recon.time_range(end)-source_recon_results.source_recon.time_range(1))/oil.enveloping.window_length);
    [~, trial_order]=sort(triallist,'ascend'); % Reorders trials into the order they were collected.    
    x=x(trial_order,:);
    xtemp=zeros(Ntrials*Nsamples_per_trial,size(x,2));
    for i=1:size(x,2)
        xtemp(:,i)= squash(repmat(x(:,i),1,Nsamples_per_trial)');
    end
    x=xtemp;
    
    pinvxtx=pinv(x'*x);
    pinvx=pinv(x);
    
    % If temporal ICA has beeen run then we want to fit temporal regressors to
    % the tICs. If spatial ICA has been run then we want to fit regressors to
    % the columns of the mixing matrix.
    
    switch oil.ica.temp_or_spat  
        case {'temporal'}
            dat=oil.ica.results.tICs(:,oil.concat_subs.results.subj_ind(sub):oil.concat_subs.results.subj_ind(sub+1)-1)';
        case {'spatial'}
            dat=oil.ica.results.mixing_matrix(oil.concat_subs.results.subj_ind(sub):oil.concat_subs.results.subj_ind(sub+1)-1,:);
    end

    for i = 1 : size(dat,2), % indexes over ICs
        
        if(first_level.use_robust_glm),
            [b,stats2]=robustfit(x,dat(:,i),'bisquare',4.685,'off'); % fits GLM robustly, Y=XB+E
            for c=1:length(contrast_list),
                switch lower(oil.ica_first_level.cope_type)
                case {'coape'}
                    first_level_results.cope(i,c)=contrast_list{c}'*abs(b);
                case {'cope'}
                    first_level_results.cope(i,c)=contrast_list{c}'*b;
                case {'acope'}
                    first_level_results.cope(i,c)=abs(contrast_list{c}'*b);
                otherwise
                    error('Invalid cope_type');                    
                end
                first_level_results.stdcope(i,c)=sqrt(contrast_list{c}'*stats2.covb*contrast_list{c});
            end;
            pe(t,:)=b; % num_timepoints x  num_pes
            
        else
            [copeout, varcopeout, coapeout, dof,peout]=glm_fast_for_meg(dat(:,i),x,pinvxtx,pinvx,contrast_list,0);
            
            switch lower(oil.ica_first_level.cope_type)
                case {'coape'}
                    first_level_results.cope(i,:)=coapeout;    
                case {'cope'}
                    first_level_results.cope(i,:)=copeout;                
                case {'acope'}
                    first_level_results.cope(i,:)=abs(copeout);
                otherwise
                    error('Invalid cope_type');                    
            end
            first_level_results.stdcope(i,:)=sqrt(varcopeout);
            pe(i,:)=peout; % num_timepoints x  num_pes
        end;
        
    end;
    
    % store single subject results in the OAT structure
    
    first_level_results.x=x;
    first_level_results.contrasts=first_level.contrast;
    first_level_results.S=first_level;
    first_level_results.source_recon=source_recon_results.source_recon;
    first_level_results.level=1;
    first_level_results.name=first_level.name;
    first_level_results.recon_method=source_recon_results.recon_method;
    first_level_results.subject_name=source_recon_results.session_name;
    
    oil.ica_first_level.results{sub}=first_level_results;
 
end
end

