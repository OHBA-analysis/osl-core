function [gamma_epoched, Ds]=teh_epoch_statetcs(hmm,D_epoched_files)

is_concat_epoched = false;

% need to make sure hmm.data_files and opt.results.spm_files_epoched line
% up

for subnum=1:length(hmm.data_files)
    D=spm_eeg_load(hmm.data_files{subnum});
    D_epoched=spm_eeg_load(D_epoched_files{subnum});
    
    sub_statepath = hmm.statepath(hmm.subj_inds==subnum);
    sub_gamma = hmm.gamma(hmm.subj_inds==subnum,:);
                
    NK=size(sub_gamma,2);
    
    % select only good samples 
    good_samples = ~all(badsamples(D,':',':',':'));
    good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);   
    
    % insert into D (continuous) object
    hmm_class=zeros(1,size(D,2),size(D,3));
    hmm_class(1,find(good_samples),1)=reshape(sub_statepath,[1,numel(find(good_samples)),1]);
    
    hmm_class_probs=zeros(NK,size(D,2),size(D,3));
    hmm_class_probs(:,find(good_samples),1)=reshape(sub_gamma',[NK,numel(find(good_samples)),1]);
 
    % add channel
    Sc=[];
    Sc.D=D;
    Sc.newchandata=[hmm_class; hmm_class_probs];
    Sc.newchanlabels{1}='Class';
    Sc.newchantype{1}='CLASS';

    for kk=1:NK,
        Sc.newchanlabels{kk+1}=['ClassPr' num2str(kk)];
        Sc.newchantype{kk+1}='CLASSPR';
    end;

    [ Dnew ] = osl_concat_spm_eeg_chans( Sc );
    
    %%%%%%%%%%
    disp('Epoching...');

    try,
        S2=D_epoched.epochinfo;
    catch
        error('No epoch info available');
    end;
    S2.D = Dnew;
    S2.epochinfo.padding = 0;
    S2.save=0;
    S2.reviewtrials=0;
    S2.bc=0; 
    disp('Doing no within-trial baseline correction at the point of epoching');

    Dnew_epoched = spm_eeg_epochs(S2);   
    
    %% get bad channels and trials from passed in D_epoched    
    if(~isempty(D_epoched) && ~isfield(Dnew_epoched, 'parcellation')),
        rejtmp=badtrials(D_epoched);
        rej=zeros(1,D_epoched.ntrials);
        rej(rejtmp)=1;
        rej(D_epoched.badtrials)=1;
        Dnew_epoched = badtrials(Dnew_epoched, 1:length(rej), rej);  

        if(length(D_epoched.badchannels)>0),
            Dnew_epoched = badchannels(Dnew_epoched, D_epoched.badchannels, ones(length(D_epoched.badchannels),1));
        end;
        Dnew_epoched.save;
    end;
    
    Ds{subnum}=Dnew_epoched;
    
    gamma_epoched(subnum,:,:)=mean(Dnew_epoched(find(strcmp(Dnew.chantype,'CLASSPR')),:,:),3);
end