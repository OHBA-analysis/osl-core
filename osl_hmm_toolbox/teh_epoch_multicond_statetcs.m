function [gamma_epoched, Ds]=teh_epoch_multicond_statetcs(hmm,epoched_cond_input_files, cont_cond_input_files)

% need to make sure hmm.data_files and opt.results.spm_files_epoched line
% up

clear gamma_epoched;

for subnum=1:length(hmm.data_files)
    
    sub_statepath = hmm.statepath(hmm.subj_inds==subnum);
    sub_gamma = hmm.gamma(hmm.subj_inds==subnum,:);
                
    NK=size(hmm.gamma,2);
    sub_statepath_cond=cell(size(epoched_cond_input_files,2),1);
    
    to=0;
    for cc=1:size(epoched_cond_input_files,2)   
       D_epoched_cond=spm_eeg_load(epoched_cond_input_files{subnum,cc});  
       from=to+1;
       to=from+size(D_epoched_cond,2)*size(D_epoched_cond,3)-1;
       
       sub_statepath_cond{cc}=sub_statepath(from:to);
       sub_gamma_cond{cc}=sub_gamma(from:to,:);
    
       D_cont_cond=spm_eeg_load(cont_cond_input_files{subnum,cc});  

        % select only good samples 
        goodsamples = good_samples(D_cont_cond);
        goodsamples = reshape(goodsamples,1,D_cont_cond.nsamples*D_cont_cond.ntrials);   

        % insert into D (continuous) object
        hmm_class=zeros(1,size(D_cont_cond,2),size(D_cont_cond,3));
        hmm_class(1,find(goodsamples),1)=permute(sub_statepath_cond{cc},[2,1,3]);
        hmm_class_probs=zeros(NK,size(D_cont_cond,2),size(D_cont_cond,3));
        hmm_class_probs(:,find(goodsamples),1)=permute(sub_gamma_cond{cc},[2,1,3]);

        % add channel
        Sc=[];
        Sc.D=D_cont_cond;
        Sc.newchandata=[hmm_class; hmm_class_probs];
        Sc.newchanlabels{1}='Class';
        Sc.newchantype{1}='CLASS';

        for kk=1:NK,
            Sc.newchanlabels{kk+1}=['ClassPr' num2str(kk)];
            Sc.newchantype{kk+1}='CLASSPR';
        end;

        Dnew = osl_concat_spm_eeg_chans( Sc );

        %%%%%%%%%%
        disp('Epoching...');
        Dtmp=spm_eeg_load(epoched_cond_input_files{subnum,cc});
        
        S2=Dtmp.epochinfo;        
        S2.D = Dnew;
        %S2.epochinfo.padding = 0;
        S2.save=0;
        S2.reviewtrials=0;
        S2.bc=0; 
        disp('Doing no within-trial baseline correction at the point of epoching');

        Dnew_epoched = spm_eeg_epochs_osl(S2);   

        %% get bad channels and trials from passed in D_epoched    
        if(~isempty(D_epoched_cond) && ~isfield(Dnew_epoched, 'parcellation')),
            rejtmp=badtrials(D_epoched_cond);
            rej=zeros(1,D_epoched_cond.ntrials);
            rej(rejtmp)=1;
            rej(D_epoched_cond.badtrials)=1;
            Dnew_epoched = badtrials(Dnew_epoched, 1:length(rej), rej);  

            if(length(D_epoched_cond.badchannels)>0),
                Dnew_epoched = badchannels(Dnew_epoched, D_epoched_cond.badchannels, ones(length(D_epoched_cond.badchannels),1));
            end;
            Dnew_epoched.save;
        end;

        Ds{subnum,cc}=Dnew_epoched;

        gamma_epoched{cc}(subnum,:,:,:)=Dnew_epoched(find(strcmp(Dnew_epoched.chantype,'CLASSPR')),:,:);
    end
end