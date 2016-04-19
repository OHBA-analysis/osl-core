function [Ds, flips, scorepath] = osl_resolve_sign_ambiguity(S)

% [D, flips, scorepath] = osl_resolve_sign_ambiguity(S)
%
% Calls findflip on passed in SPM objects and applies the flips to resolve sign ambiguity
%
% INPUTS:
%
% S.Ds          - cell array of filenames (with full paths) to 
%               SPM objects of PARCELATED source space MEEG object
%               NOTE that this function assumes that all parcels are good
%               (i.e. it does not check for bad parcels using
%               D.badchannels) and that parcels have  'VE' as a chantype.
%
% S.options     - specifies options to pass to the func findflip. 
%               Default: uses 
%               findflip defaults.
%
% S.prefix      - prefix for new sign flipped SPM files. Default: 'sf_'
%
% OUTPUTS:
%
% Ds            - new MEEG object containing the flipped time courses
% flips         - (num_Ds x num_parcels binary matrix indicating the sign
%               flipping that has been done.
% scorepath     - output from internal call to findflip
%
% MWW 2016

global OSLDIR 

num_subj=length(S.Ds);

if isfield(S,'options')
    options = S.options;
else
    options=[]; 
end;

if ~isfield(S,'prefix')
    S.prefix = 'sf_'
end

%% extract data and concatenate over subjects

subj_indexes = [];
parceldata = [];
Ts = [];
for subj=1:num_subj,
    
    D = spm_eeg_load(S.Ds{subj});

    if subj==1
        chanind = strmatch('VE', D.chantype); % 'VE' is chantype for parcels
    else
        chanind_new = strmatch('VE', D.chantype);
        if length(chanind) ~= length(chanind_new)
            error('All S.Ds must correspond to the same parcellation');
        end
    end
    
    subjparceldata = reshape(D(chanind,:,:),[length(chanind),D.nsamples*D.ntrials]);
    good_samples = ~all(badsamples(D,':',':',':'));
    good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
    subjparceldata = subjparceldata(:,good_samples);        
    
    parceldata = [parceldata subjparceldata];
    subj_indexes = [subj_indexes subj*ones(1,size(subjparceldata,2))]; 
    Ts=[Ts size(subjparceldata,2)];
end

%%%%%%%%
%% do the actual work

% parceldata': a (no. of time points X no. of channels) matrix containing the time series
% Ts: a vector containing the length of each trial/segment/subject (the sum of T must be equal to the number of rows of X).

[flips,scorepath] = findflip(parceldata',Ts,options);

%flipped_parceldata = flipdata(parceldata',Ts,flips)';

%%%%%%%%
%% apply flips to copy of spm objs

clear Ds;
for subj=1:num_subj,
    
    S2=[];
    S2.D=S.Ds{subj};
    [pathstr,filestr] = fileparts(S2.D);
    S2.outfile=[pathstr '/' S.prefix filestr];
    Dnew = spm_eeg_copy(S2);

    % Apply flips
    apply_flips=-flips(subj,:)*2+1;
    Dnew(chanind,:,:) = Dnew(chanind,:,:).*repmat(apply_flips',[1,size(Dnew,2)]);
    Dnew.save;
    
    Ds{subj}=Dnew;
    
end

end