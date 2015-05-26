function fif_out = osl_call_maxfilter(S)
% Calls Elekta MaxFilter to apply signal space separation (SSS) to input
% .fif file.
%
% Syntax: fif_out = osl_call_maxfilter(S)
%
% S needs to contain:
%   -   fif: full name of input fif file (inc. path if not in present
%       working directory, with or without .fif file extention.
%
%   -   fif_out (optional): manually specify name of output file.
%
%   -   logfile (optional): set to 1 to produce log file of SSS output. Log 
%       file name is fif_out_log.txt. Set to 0 for no log file (default).
%
%   -   downsample_factor (optional): integer number >1. Suggested value of 
%       4 for 1000KHz data. If S.downsample_factor is not set then no 
%       downsampling occurs.
%
%   -   spmfile (optional): Include the name of corresponding SPM format 
%       data. Maxfilter will use bad channel information contained in SPM 
%       object to ignore bad channels.
%
%   -   nosss (optional): set to 1 to call Maxfilter without applying SSS.
%
%   -   maxfilt_dir (recommended): explicity tell OMT where to find
%       MaxFilter. Defaults to S.maxfilt_dir = '/neuro/bin/util'.
%
%   -   st (optional): set to st.do to 1 to use temporal_extension or 0 
%       (default) to not use it. Must be combined with SSS.
%
%   -   movement_compensation (optional): set to 1 to use cHPI movement
%       compensation or 0 (default) to not use it. Must be combined with SSS.
%
%   -   movecomp_call (optional, default = ''): string specifying a manual 
%       movement compensation call. Setting this field overrides 
%       movement_compensation. Must be combined with SSS.
%
%   -   autobad_off (optional): set to 1 to switch off MaxFilter's automated
%       badchannel detection. Default = 0.
%
%   -   trans_ref_file (optional, default = ''): allow transformation to the 
%       specified reference .fif file using the -trans maxfilter option
%
%   -   headpos (optional, default = 0): set to 1 to output head position 
%       information from the recording to a text file.
%
%   -   bad_epochs (optional): used in -skip call to maxfilter. Should be a
%       matrix of dimensions num_bad_epochs x 2, where each row is the start 
%       and stop time of the epoch to skip.'), outputs
%
%   -   cal_file (optional):  full path to a fine calibration file to use
%
%   -   ctc_file (optional):  full path to a cross-talk matrix file to use
%
% OUTPUTS:
%   - fif_out: the name of the output file.
%
% HL 061011


if ~isfield(S, 'fif')
    error('No .fif file specified.')
else
    [pathstr,filestr,~] = fileparts(S.fif);
    S.fif = fullfile(pathstr,filestr);
end

if ~isfield(S,'fif_out')
    fif_out=[S.fif '_sss'];
else
    fif_out=S.fif_out;
end

if isfield(S,'downsample_factor')
    if S.downsample_factor==1
        ds_call='';
    else
        ds_call=[' -ds ' num2str(S.downsample_factor)];
    end
else
    ds_call='';
end

if isfield(S,'spmfile')
    D=spm_eeg_load(S.spmfile);
    badchans=[]; chanind=D.badchannels;
    for i=1:length(D.badchannels)
        chanlab=D.chanlabels(chanind(i));
        badchans=[badchans ' ' chanlab{1}(4:end)];
    end
    if ~isempty(badchans)
    badchan_call = [' -bad ' badchans];
    else
        badchan_call = '';
    end
else
    badchan_call='';
end

if isfield(S,'cal_file')
    cal_call = [' -cal ' S.cal_file ];
else
    cal_call = '';
end

if isfield(S,'ctc_file')
    ctc_call = [' -ctc ' S.ctc_file ];
else
    cal_call = '';
end

if isfield(S,'nosss')
    if S.nosss==1
        nosss_call=[' -nosss'];
        if ~isfield(S,'fif_out')
            fif_out=[S.fif '_nosss'];
        end
    else
        nosss_call='';
    end
else
    nosss_call='';
end
if isfield(S,'maxfilt_dir')
    max_dir=[S.maxfilt_dir];
else
    max_dir='/neuro/bin/util';
end
disp(['Assuming Maxfilter found in ' max_dir]);

if isfield(S,'logfile')
    if S.logfile==1
    log_call=[' >& ' fif_out '_log.txt'];
    else
      log_call='';
    end
else
    log_call='';
end
if isfield(S,'movecomp_call') % ICG ADDITION, ALLOW CUSTOM MOVECOMP CALLS
    movecomp_call=S.movecomp_call;
    if isfield(S,'movement_compensation') && S.movement_compensation==0
        warning('Movement compensation call is specifed, even though movemen_compensation == 0!');
    end
elseif isfield(S,'movement_compensation')
    if S.movement_compensation==1 && isempty(nosss_call)
        movecomp_call=[' -movecomp'];
    else
        movecomp_call='';
    end
else
    movecomp_call='';
end
if isfield(S,'trans_ref_file') % ICG ADDITION, ALLOW TRANFORM TO REFERENCE FILE
    trans_call= [' -trans ' S.trans_ref_file '.fif'];
else
    trans_call='';
end
if isfield(S,'headpos') % ICG ADDITION, ALLOW OUTPUT OF HEAD MOVEMENT INFO
    if S.headpos==1
        headpos_call = [' -hp ' fif_out '_headpos.txt'];
    else
        headpos_call='';
    end
else 
    headpos_call='';
end
if isfield(S,'autobad_off')
    if S.autobad_off==1 && isempty(nosss_call)
        autobad_call=[' -autobad off'];
    else
        autobad_call='';
    end
else
    autobad_call='';
end

st_call='';
if isfield(S,'st') % SSS with ST
    if S.st.do % SSS with ST
        if ~isfield(S.st,'win') | isempty(S.st.win),
            %error('specify S.st.win') % use 4 as default?
            S.st.win=4;
        end;
        
        if ~isfield(S.st,'corr') | isempty(S.st.corr),
            %error('specify S.st.corr') % use .98 as default?
            S.st.corr=0.98;            
        end
        st_call= sprintf('  -st %d -corr %g',S.st.win,S.st.corr);
    end;
end

skip_call='';
if isfield(S,'bad_epochs'), % -skip option
    if ~isempty(S.bad_epochs),        
        skip_call= [' -skip'];

        if size(S.bad_epochs,2)~=2,
            error('S.bad_epochs should be num_bad_epochs x 2, where each row is the start and stop time of the epoch to skip.');
        end;
        
        for ee=1:size(S.bad_epochs,1),
            BadEpoch=S.bad_epochs(ee,:);

            %if BadEpoch(1)==-1, BadEpoch(1)=D.time(1); end;        
            %if BadEpoch(2)==-1, BadEpoch(2)=D.time(end); end;
            skip_call= [skip_call ' ' num2str(BadEpoch(1)) ' ' num2str(BadEpoch(2))];
        end;   
    
    end; 
end;

rmcall=['rm -f ' fif_out];
disp(['Deleting file '  fif_out]);
runcmd(rmcall);
rmcall=['rm -f ' fif_out '.fif'];
runcmd(rmcall);

maxfilter_call=[max_dir '/maxfilter -f ' S.fif '.fif -o ' fif_out '.fif'  nosss_call ds_call badchan_call movecomp_call st_call autobad_call trans_call headpos_call skip_call ' -format float -v ' cal_call ctc_call log_call]

try,
    runcmd(maxfilter_call)
catch
    
end;
