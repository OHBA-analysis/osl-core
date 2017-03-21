function [D fname fig_handles fig_names] = osl_convert_script(Sin)

% [D fname] = osl_convert_script(Sin)
% 
% Converts .fif file to continuous spm data file
%
% Mandatory Inputs:
%
% Sin.fif_file : file name of input fif file
% Sin.spm_file : file name of outputted spm file
%                
% Optional inputs:
%
% Sin.num_condition_trials_expected : number of expected trials for each condition 
% Sin.condition_code : list of condition trigger codes 
%
% MW

try; tmp=Sin.fif_file; catch, error('fif_file not specified');end;
try; tmp=Sin.spm_file; catch, error('spm_file not specified');end;
try; tmp=Sin.trigger_channel_mask; catch, Sin.trigger_channel_mask='0000000001111111';end;
%Sin.trigger_channel_mask='0000000000000111';

fig_handles=[];
fig_names={};

S = [];
S.dataset = Sin.fif_file;
S.outfile = Sin.spm_file;


% Create directory if it doesn't exist
pathstr = fileparts(Sin.spm_file);
if ~isdir(pathstr)
    mkdir(pathstr);
end


[dir nam ext]=fileparts(Sin.fif_file);

if strcmp(ext,'.ds')
    %S.channels=ctf_chans;     % HL Mod 1.2 - CTF channel selection
    %TC 2013 the above works only for a 275 channel ctf system
    %To make more general, we want to read stim channels, the clock channel
    %and meg grads. Ideally data should actually be beamformed in 3rd gradient
    %which would require the reference channels but that would require
    %modification of the forward model (and you would also need to read in
    %balancing coefs). For now the data is read in unbalanced (raw) - this
    %does not take advantage of the noise corrections available to CTF
    %systems.
    S.channels='ctf_grads';    
    S.checkboundary = 0;  
    
    if isfield(Sin, 'artefact_channels'),
        S.artefact_channels = Sin.artefact_channels;
    else
        S.artefact_channels = [];
    end
    isctf=1;
else
    S.channels = 'all';
    S.checkboundary = 1;
    S.fixoxfordneuromag=Sin.trigger_channel_mask;
    isctf=0;
end

convertfun = @spm_eeg_convert_4osl; 

S.usetrials = 1;
S.datatype = 'float32-le';
S.eventpadding = 0;
S.saveorigheader = 1;
S.conditionlabel = {'Undefined'};
S.inputformat = [];
S.mode='continuous';

D = feval(convertfun,S);


%% Add coil to channel mapping in D.sensors (used later for local spheres)
if ~isempty(D.sensors('MEG'))
    sens = D.sensors('MEG');
    sens.coilchan = sens.tra==1;
    D = sensors(D,'MEG',sens);
end

%%%%%%%%%%%%%%%%%%%%%
%% MWW added back in from osl1.2.beta.15 - Oct 2013
if 1,
    eve=events(D,1);    

    if(length(eve)>0)

        setval=0;

        for i=1:length(eve),

            if(eve(i).value>0)
                vals(i)=eve(i).value;  
            else
                eve(i).value=1;
                vals(i)=eve(i).value;
                setval=1;
            end;
        end;

        if(setval)
            D=events(D,1,eve);
            D.save;
            warning('Some events detected with no value available in D.events, so setting them to 1');
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%

if ~isctf
    eve=events(D,1);    

    if(length(eve)>0)

        for i=1:length(eve),                  
            vals(i)=eve(i).value;   
        end;

        bins=unique(vals); % HL mod 1.1
        bins=[bins(1)-1 bins bins(end)+1]; % HL mod 1.1
        h=hist(vals,bins); % HL mod 1.1
        fig_handles(1)=sfigure;bar(bins,h/2); % HL mod 1.1, correct by 2 due to up and down triggers
        fig_names{1}='events_hist';
        title('Histogram of events corrected for button presses');xlabel('Trigger codes');ylabel('# of triggers');

        if(isfield(Sin,'num_condition_trials_expected'))
            if(~isempty(Sin.num_condition_trials_expected))

                % check to see if correct number of events have been found
                uvals=unique(vals);

                    wrong=0;
                for i=1:length(Sin.condition_code),
                    ind=find(uvals==Sin.condition_code(i));

                    xs(i)=ind;
                    ys(i)=h(ind);

                    if(ys(i)~=Sin.num_condition_trials_expected(i))
                        wrong=1;
                    end;
                end;

                sfigure;bar(1:4, ys-Sin.num_condition_trials_expected);

                if(wrong)
                   warning('Wrong number of events');
                   print('-dpng',[working_dir '/' S.outfile '_num_trials']);

                   h=hist(vals,unique(vals));
                   sfigure;bar(unique(vals),h);

                   print('-dpng',[working_dir '/' S.outfile '_num_trials2']);

                else
                   disp('Correct number of events');
                end;
            end;
        end;

    else
        warning('No events detected');
    end;

    D=events(D,1,eve);
end

D.save;
fname=[D.path '/' D.fname];

end
