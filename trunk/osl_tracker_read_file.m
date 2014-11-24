function con = osl_tracker_read_file(S)
% TRACKER_READ_FILE read in the tab delimitated report from Eye Link Viewer
%
% Input:
%   files (required)    : Tracker sample report:
%                         Column headers: TIMESTAMP, LEFT_GAZE_X, LEFT_GAZE_Y, LEFT_IN_BLINK,
%                         LEFT_IN_SACCADE, RIGHT_GAZE_X, RIGHT_GAZE_Y, RIGHT_IN_BLINK, RIGHT_IN_SACCADE
%   trunc (optional)    : [start, end] Remove samples before start time (s) and after end time
%   blink_win (required): Blink Window [before, after] Samples before and after blink to ignore/avg
%   D (alternative to trunk) : SPM object with eyetracker channel for
%                              estimating lags.
%   eyechan (needed with D)  : User specified channel in D that has
%                              eyetracker recording.
%
% Output:
%   out: con.t, con.l_x, con.l_y, con.l_blink, con.l_blink_edge, con.l_blink_dur
%        con.l_sac, con.l_sac_edge, con.l_sac_dur, con.l_sac_r, con.l_sac_theta
%        con.r_x, con.r_y, con.r_blink, con.r_blink_edge, con.r_blink_dur
%        con.r_sac, con.r_sac_edge, con.r_sac_dur, con.r_sac_r, con.r_sac_theta
%
% Examples:
%   continuous_position = tracker_read_file(file,[244.297,12.445], [100,100]);
%   continuous_position = tracker_read_file(S);
%   continuous_position = tracker_read_file(file,[],[100,100],D,eyechan);
% TODO: Expand to binocular to use average left/right gaze and compare
% left/right saccade/blink. This will help reduce false blinks when the
% tracker looses an eye. 
%
% Tyler Ferro (University of oxford, MSc Biomeical Engineering, OHBA, 2011)
% 
%
%
%
% Changelog: 
% 
% AB fix 1 - ensure consistency when removing last sample
% AB fix 2 - normalisation of xcorr output
% AB fix 3 - remove hard coding of eyetracker sampling rate
%


%% Accept single structure input

if isfield(S,'file')
    file=S.file;
else
    error('No eyetracker file specified')
end
if isfield(S,'blink_win')
    blink_win=S.blink_win;
else
    blink_win=[100 100];
end    
if isfield(S,'trunc')
    trunc=S.trunc;
else
    trunc=[];
end  
if isfield(S,'D')
    D=S.D;
end
if isfield(S,'eyechan')
    eyechan=S.eyechan;
end  

if (~exist('eyechan','var') && ~exist('D','var')) && ~exist('trunc','var')
    if (exist('eyechan','var') && ~exist('D','var'))
        error('SPM object not provided - Insufficient information provided to align eyetracker to MEG data') 
    elseif (~exist('eyechan','var') && exist('D','var'))
        error('Eyetracker channel in SPM object not specified - Insufficient information provided to align eyetracker to MEG data')
    else
        error('Insufficient information provided to align eyetracker to MEG data')
    end
end
%% Read File

fprintf('\nProcessing tracker file\n ')

% Initialize structure 
% init = zeros(1,length(lines-1));
con = struct('t',[],'l_x',[],'l_y',[],'l_blink',[],'l_sac',[],'r_x',[],'r_y',[],'r_blink',[],'r_sac',[]);

% Parse file
fid = fopen(file,'rt');
hdr = fgetl(fid);
column_headers=textscan(hdr,'%s%s%s%s%s%s%s%s%s','Delimiter','\t');
str=[];
for i=1:length(column_headers)
    str=[str ', ' cell2mat(column_headers{i})];
end
fprintf('\n%s%s%s\n','Detected columns with the following eyetracker variables: ',str,'.')

clear lin;

iter = 1;
lines = 1;
               
while ~feof(fid)
    lin = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'Delimiter','\t','TreatAsEmpty', '.');   
    for col_ind=1:length(column_headers)
        variable=cell2mat(lower(column_headers{col_ind}));
        if ~isempty(variable)
            switch variable
                case 'timestamp'
                con.t        =   [con.t; lin{col_ind}];
                case 'left_gaze_x'
                con.l_x      =   [con.l_x; lin{col_ind}];
                case 'left_gaze_y'
                con.l_y      =   [con.l_y; lin{col_ind}];
                case 'left_in_blink'
                con.l_blink  =   [con.l_blink; lin{col_ind}];
                case 'left_in_saccade'
                con.l_sac    =   [con.l_sac; lin{col_ind}];
                case 'right_gaze_x'
                con.r_x      =   [con.r_x; lin{col_ind}];
                case 'right_gaze_y'
                con.r_y      =   [con.r_y; lin{col_ind}];
                case 'right_in_blink'
                con.r_blink  =   [con.r_blink; lin{col_ind}];
                case 'right_in_saccade'
                con.r_sac    =   [con.r_sac; lin{col_ind}];
            end
        end
    end
    iter = iter + 1;
    lines = lines + length(lin{1});
end

fclose(fid);

%% Adjust time scale

% Remove last sample, sometimes last timestamp is NaN
con = structfun(@(x) x(1:end-1),con,'UniformOutput',false); % AB fix 1

% Eyetracker sample frequency - AB fix 3
fs_con = 1000/diff(con.t(1:2));

% Convert timescale from milliseconds to seconds
con.t = con.t*10^-3;

% Start time to zero
con.t = con.t - con.t(1);




%% Estimate Lags

if isempty(trunc) % Estimate lags automatically
    
    fprintf('%s\n','Automatically estimating lags from user-provided SPM object')
    fif_track=interp(D(eyechan,:),fs_con/D.fsample);

    if ~isempty(con.l_blink)
        [l_C,l_LAGS]=xcorr(con.l_blink,fif_track);
        l_C = l_C /(norm(con.l_blink)*norm(fif_track)); % AB fix 2
        
        [l_m l_ind]=max(abs(l_C));
    end
    if ~isempty(con.r_blink)
        [r_C,r_LAGS]=xcorr(con.r_blink,fif_track);
        r_C = r_C /(norm(con.r_blink)*norm(fif_track)); % AB fix 2
        [r_m r_ind]=max(abs(r_C));
    end
    
    if r_m>l_m
        delay=r_LAGS(r_ind);
        et_blink=con.r_blink;
        fprintf('%s\n%s\n','Eyetracker channel in FIF best matches right eye.', 'Using right eye blinks for lag estimation.')
    else
        delay=l_LAGS(l_ind);
        et_blink=con.l_blink;
        fprintf('%s\n%s\n','Eyetracker channel in FIF best matches left eye.', 'Using left eye blinks for lag estimation.')
    end
    
    if delay<0
        fprintf('%s','Delay is negative; Eye-tracker started recording after fif. Padding eyetracker data accordingly.')
        padding=zeros(abs(delay),1);
        con.t   =  [(delay:1:-1)'/fs_con; con.t];
        con.t = con.t - con.t(1);
        con.l_x   =  [padding; con.l_x];
        con.l_y   =   [padding;con.l_y];
        con.l_blink = [padding;con.l_blink];
        con.l_sac =   [padding;con.l_sac];
        con.r_x   =   [padding;con.r_x];
        con.r_y   =   [padding;con.r_y];
        con.r_blink = [padding;con.r_blink];
        con.r_sac =   [padding;con.r_sac];
        delay=1;
    else
    fprintf('%s%0.5g%s\n','Eyetracker channel is ', con.t(delay),'s ahead of the fif recording.')
    end
    trunc_index=[delay delay+length(fif_track)-1];
    
else % User manually provides lags
    % Find index of times to remove
    %       uint32(*10000) because matlab doesn't like to find decimals
    trunc_index(1) = find(uint32(con.t*10000) == trunc(1)*10000);
    trunc_index(2) = find(uint32(con.t*10000) == (uint32(con.t(end)) - trunc(2))*10000);
end

%% Remove samples
if trunc_index(2)>length(con.t)
    fprintf('%s',' Eye-tracker stopped recording before fif. Padding eyetracker data accordingly.')
        con.t   = [con.t; ((length(con.t)+1):trunc_index(2))'/1000];
        con.l_x(end+1:trunc_index(2))   =  0;
        con.l_y(end+1:trunc_index(2))   =  0;
        con.l_blink(end+1:trunc_index(2))   =  0;
        con.l_sac(end+1:trunc_index(2))   =  0;
        con.r_x(end+1:trunc_index(2))   =  0;
        con.r_y(end+1:trunc_index(2))   =  0;
        con.r_blink(end+1:trunc_index(2))   =  0;
        con.r_sac(end+1:trunc_index(2))   =  0;
        trunc_index(2)=length(con.t);
end

con.t   =   con.t(trunc_index(1):trunc_index(2));
con.l_x   =   con.l_x(trunc_index(1):trunc_index(2));
con.l_y   =   con.l_y(trunc_index(1):trunc_index(2));
con.l_blink = con.l_blink(trunc_index(1):trunc_index(2));
con.l_sac =   con.l_sac(trunc_index(1):trunc_index(2));
con.r_x   =   con.r_x(trunc_index(1):trunc_index(2));
con.r_y   =   con.r_y(trunc_index(1):trunc_index(2));
con.r_blink = con.r_blink(trunc_index(1):trunc_index(2));
con.r_sac =   con.r_sac(trunc_index(1):trunc_index(2));

% Start time to zero
con.t = con.t - con.t(1);

%figure; plot(D.time,D(eyechan,:));
%ho; plot(con.t,con.r_blink,'r--');
%ho; plot(con.t,con.l_blink,'g-.');
%xlabel('Time (s)'); title('Fif, right and left blink time courses after alignment'); legend('Fif', 'Eyetracker - right','Eyetracker - left');

%% Init

con.l_blink_dur = zeros(length(con.t),1);
con.l_sac_dur = zeros(length(con.t),1);
con.l_sac_theta = zeros(length(con.t),1);
con.l_sac_r = zeros(length(con.t),1);
con.r_blink_dur = zeros(length(con.t),1);
con.r_sac_dur = zeros(length(con.t),1);
con.r_sac_theta = zeros(length(con.t),1);
con.r_sac_r = zeros(length(con.t),1);
%% Extract Blink
if ~isempty(con.l_blink)
con.l_blink(1) = 0; con.l_blink(end) = 0;
con.l_blink_edge = [0; diff(con.l_blink)];
con.l_blink_dur(con.l_blink_edge > 0) = con.t(con.l_blink_edge < 0) - con.t(con.l_blink_edge > 0);
end
if ~isempty(con.r_blink)
con.r_blink(1) = 0; con.r_blink(end) = 0;
con.r_blink_edge = [0; diff(con.r_blink)]; 
con.r_blink_dur(con.r_blink_edge > 0) = con.t(con.r_blink_edge < 0) - con.t(con.r_blink_edge > 0);
end
%% Replace blink x/y with avg before/after

% Left Eye
if ~isempty(con.l_blink)
l_neg_edge = find(con.l_blink_edge < 0);ind_neg=find(l_neg_edge>blink_win(1)&l_neg_edge<length(con.l_x)-blink_win(2));
l_pos_edge = find(con.l_blink_edge > 0);ind_pos=find(l_pos_edge>blink_win(1)&l_pos_edge<length(con.r_x)-blink_win(2));
l_neg_edge=l_neg_edge(intersect(ind_neg,ind_pos));
l_pos_edge=l_pos_edge(intersect(ind_neg,ind_pos));
for iter = 1:length(l_pos_edge)
    con.l_x(l_pos_edge(iter)-blink_win(1):l_neg_edge(iter)+blink_win(2)) = (con.l_x(l_pos_edge(iter)-blink_win(1)) + con.l_x(l_neg_edge(iter)+blink_win(2)))/2;
    con.l_y(l_pos_edge(iter)-blink_win(1):l_neg_edge(iter)+blink_win(2)) = (con.l_y(l_pos_edge(iter)-blink_win(1)) + con.l_y(l_neg_edge(iter)+blink_win(2)))/2;
    con.l_sac(l_pos_edge(iter)-blink_win(1):l_neg_edge(iter)+blink_win(2)) = 0;
end
end
% Right Eye
if ~isempty(con.r_blink)
r_neg_edge = find(con.r_blink_edge < 0);ind_neg=find(r_neg_edge>blink_win(1)&r_neg_edge<length(con.l_x)-blink_win(2));
r_pos_edge = find(con.r_blink_edge > 0);ind_pos=find(r_pos_edge>blink_win(1)&r_pos_edge<length(con.r_x)-blink_win(2));
r_neg_edge=r_neg_edge(intersect(ind_neg,ind_pos));
r_pos_edge=r_pos_edge(intersect(ind_neg,ind_pos));
for iter = 1:length(r_pos_edge)
    con.r_x(r_pos_edge(iter)-blink_win(1):r_neg_edge(iter)+blink_win(2)) = (con.r_x(r_pos_edge(iter)-blink_win(1)) + con.r_x(r_neg_edge(iter)+blink_win(2)))/2;
    con.r_y(r_pos_edge(iter)-blink_win(1):r_neg_edge(iter)+blink_win(2)) = (con.r_y(r_pos_edge(iter)-blink_win(1)) + con.r_y(r_neg_edge(iter)+blink_win(2)))/2;
    con.r_sac(r_pos_edge(iter)-blink_win(1):r_neg_edge(iter)+blink_win(2)) = 0;
end
end
%% Extract Saccade
% Left Eye
if ~isempty(con.l_sac)
con.l_sac(1) = 0; con.l_sac(end) = 0;
con.l_sac_edge = [0; diff(con.l_sac)];
con.sac_dur(con.l_sac_edge > 0) = con.t(con.l_sac_edge < 0) - con.t(con.l_sac_edge > 0);

l_dx = con.l_x(con.l_sac_edge < 0) - con.l_x(con.l_sac_edge > 0);
l_dy = con.l_y(con.l_sac_edge < 0) - con.l_y(con.l_sac_edge > 0);

[l_theta,l_r] = cart2pol(l_dx,l_dy);
con.l_sac_theta(con.l_sac_edge > 0) = l_theta;
con.l_sac_r(con.l_sac_edge > 0) = l_r;
end

% Right Eye
if ~isempty(con.r_sac)
con.r_sac(1) = 0; con.r_sac(end) = 0;
con.r_sac_edge = [0; diff(con.r_sac)];
con.sac_dur(con.r_sac_edge > 0) = con.t(con.r_sac_edge < 0) - con.t(con.r_sac_edge > 0);

r_dx = con.r_x(con.r_sac_edge < 0) - con.r_x(con.r_sac_edge > 0);
r_dy = con.r_y(con.r_sac_edge < 0) - con.r_y(con.r_sac_edge > 0);

[r_theta,r_r] = cart2pol(r_dx,r_dy);
con.r_sac_theta(con.r_sac_edge > 0) = r_theta;
con.r_sac_r(con.r_sac_edge > 0) = r_r;
end
fprintf('- Done\n');

end

