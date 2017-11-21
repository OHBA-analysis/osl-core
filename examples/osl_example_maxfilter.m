%% HOW TO MAXFILTER

%%
% |maxfilter| is a program provided by Elekta, which implements a spatial signal
% space separation (SSS) algorithm to remove external noise from MEG
% recordings. To oversimplify, it uses a beamforming-like approach to identify
% and subtract signals originating outside the scanner unit. In addition to
% SSS denoising, and as part of the SSS algorithm, |maxfilter| is able to
% correct for movements and differences in head position between subjects, re-
% projecting the data onto the MEG sensors as if the data had been recorded
% with the head in a different position.

%%
% |maxfilter| is distributed and licensed by Elekta - if you have an Elekta
% scanner, then you probably already have access to |maxfilter| through your
% centre. To use |osl_maxfilter| to call |maxfilter|, you need to
% 
% * Be on a unix-like platform (not Windows) - we have not tested this on
%   Windows
% * Be running Matlab on a machine with access to the |maxfilter| binaries
% 
% If you are not at OHBA, you will need to use the optional |maxfilter_path|
% argument to specify the location of |maxfilter| on your system.
%
% |osl_maxfilter| supports three usage modes/methods. These are
%
% * |nosss| - Just read the file and perform downsampling. Downsampling is
%   important to remove HPI signals
% * |sss| - Perform standard source space separation
% * |tsss| - Perform temporal source space separation
%
% The |osl_maxfilter| function takes in the name of a '.fif' file and a
% method, and returns the name of the '.fif' file output by |maxfilter|. For
% example
fif_file = 's01_lizard.fif';
fif_out = osl_maxfilter(fif_file,'nosss'); 
fif_out = osl_maxfilter(fif_file,'sss'); 
fif_out = osl_maxfilter(fif_file,'tsss'); 

%%
% By default, the method is prepended to the filename. For example, if you pass in
fif_file = 'test.fif'
fif_out = osl_maxfilter(fif_file,'nosss'); 

%%
% the resulting output file will be called |nosss_test.fif|. You can override this by 
% specifying |fif_out|
fif_out = osl_maxfilter(fif_file,'nosss','fif_out','my_output.fif'); 

%%
% In addition to the '.fif' file, the command line output will be logged to 'my_output_log.txt'
% and 'my_output_err.txt', respectively. If you want to monitor the progress of |maxfilter|, set
% 'verbose' to |true| which will print the output of |maxfilter| to the Matlab command window 
% (as well as saving it to the log file).
fif_out = osl_maxfilter(fif_file,'nosss','verbose',true);
 
%% Double maxfilter
% The commands above are all you need to run |maxfilter|. However, if there are bad channels present, the 
% output of |maxfilter| can include artefacts. The approach usually taken is to run maxfilter, then 
% detect any bad channels, and then run maxfilter again with those channels excluded. You can do this
% by passing the bad channels to |osl_maxfilter| as follows:
%
% First, do a nosss
nosss_fif = osl_maxfilter(fif_file,'nosss');

% Then, find bad channels
D = osl_import(nosss_fif);
D = osl_detect_artefacts(D,'badtimes',false);

%%
% Now use these to maxfilter
[sss_fif,bad_segments] = osl_maxfilter(fif_file,'sss','badchannels',D.chanlabels(D.badchannels),'fif_out','awesome_sss_final.fif','verbose',true);

%%
% Note that bad times get set to zero, and these times are checked against the output log file. These
% times can be specified in |osl_import| so they are correctly treated as artefacts
D_sss = osl_import(sss_fif,'bad_segments',bad_segments);

%% Temporal SSS
% If you want to use the temporal extension to SSS, you can specify the method as 'tsss' instead of 'sss'
[sss_fif,bad_segments] = osl_maxfilter(fif_file,'tsss','badchannels',D.chanlabels(D.badchannels));
D_tsss = osl_import(sss_fif,'bad_segments',bad_segments);

