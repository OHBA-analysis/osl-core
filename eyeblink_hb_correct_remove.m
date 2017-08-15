function D=eyeblink_hb_correct_remove(S)

% eyeblink_hb_correct_remove(S)
%
% Used to examine and select spatial confounds for removal from raw M/EEG
% data
% Requires heartbeat_correact and/or eyeblink_correct to have been run
% first
% 
% subjid is e.g. meg301
% method is 'select' 'correct' or 'both' (default 'both')
% runs is list of run numbers to work on (default 1:3)
%
% LH + MW  

fname=S.D;
hbfname=S.hbfname;
blinkfname=S.blinkfname;
method=S.method;
runs=S.runs;

if strcmp(method,'select')||strcmp(method,'both')

    try
        runnum = runs;
        error
        sfigure;
        for runnum = runs
       
            svdfname = fullfile([fname '_svdoutputs'],'hb_svd.mat');
         
            load(svdfname);

            %display component magnitudes:
            comp_amp = diag(svdoutput.s); %loaded from svdfname
            subplot(3,1,runnum); hold on; 
            plot(comp_amp); 
            plot(comp_amp(1:5),'r.');
            title(['Heartbeat SVD component amplitudes - run ' num2str(runnum)]);
        end
        for runnum = runs
            disp(['Loading components for run ' num2str(runnum)]);

            %try and load heartbeat artefact, accept or reject components
        
            D = spm_eeg_load(hbfname);

            Dorig = spm_eeg_load(fname);

            S2 = [];
            S2.D = D;
            S2.Draw = Dorig;
            S2.sampwin = 1:10000;
            S2.compchan = find(strncmp(Dorig.chanlabels,'ECG',3));
            S2.method = 'reject';
            S2.fname=[S.D '_ecg'];
            S2.auto_remove = S.auto_remove;
            D = view_spatial_confounds(S2);
%            close all;
            hbcomp{runnum} = sconfounds(D);

        end
        hbok = 1;
    catch
        disp(['No heartbeat artefact files found for ' S.D ]);
        hbcomp = [];
        hbok = 0;
    end





        %try and load eyeblink artefact, accept or reject components
   % try

        %sfigure;
        %for runnum = runs
            %svdfname = [fname '_svdoutputs/blink_svd.mat'];
   
            %load(svdfname);

            %display component magnitudes:
            %comp_amp = diag(svdoutput.s); %loaded from svdfname
            
            %subplot(3,1,runnum); hold on; 
            %plot(comp_amp); 
            %plot(comp_amp(1:5),'r.');
            %title(['Eyeblink SVD component amplitudes - run ' num2str(runnum)]);
        %end

        for runnum = runs
            disp(['Loading components for run ' num2str(runnum)]);

            D = spm_eeg_load(blinkfname);

            fname = S.D;
            Dorig = spm_eeg_load(fname);

            S2 = [];
            S2.D = D;
            S2.Draw = Dorig;
            S2.sampwin = 1:10000;
            S2.method = 'reject';
            S2.compchan = find(strncmp(Dorig.chanlabels,'EOG',3));
            S2.fname=[S.D '_eb'];
            S2.auto_remove=S.auto_remove;

            D = view_spatial_confounds(S2);
            %close all;
            blinkcomp{runnum} = sconfounds(D);
        end
            blinkok = 1;
 %   catch
  %          disp(['No eyeblink artefact files found for ' S.D]);
  %          blinkcomp = [];
 %           blinkok = 0;
 %   end

    fname = fullfile([S.D '_svdoutputs'],'rejectfs.mat'); 
    
    save(fname,'blinkcomp','hbcomp','blinkok','hbok');
end


if strcmp(method,'correct')||strcmp(method,'both')
    fname = fullfile([S.D '_svdoutputs'],'rejectfs.mat');  %contains blinkcomp and hbcomp
    load(fname);
    
    for runnum = runs

        %amalgamate eyeblink and heartbeat components
        compr = [];
        if blinkok
            compr = [compr blinkcomp{runnum}];
        end    
        if hbok
            compr = [compr hbcomp{runnum}];
        end

        if ~isempty(compr)
            %remove components from original dataset
            disp('Removing components.');
            fname = S.D;
            D = spm_eeg_load(fname);

            if isempty(D.fiducials);
                error(['Need to get fiducials for file ' fname]);
            end

            ncomp = size(compr,2);
            [sel1, sel2] = spm_match_str(D.chanlabels(D.meegchannels), D.chanlabels(setdiff(D.meegchannels, D.badchannels)));
            sconf = [];
            sconf.label = D.chanlabels(D.meegchannels);
            sconf.coeff = nan(length(sconf.label), ncomp);
            sconf.coeff(sel1, :) = compr(sel2, 1:ncomp);
            sconf.bad = ones(length(sconf.label), 1);
            sconf.bad(sel1, :) = 0;
            D = sconfounds(D, sconf);

            S2 = [];
            S2.D = D;
            S2.correction = 'SSP';
              D = spm_eeg_correct_sensor_data(S2);
            D.save; 
        else
  
            
            disp('No components rejected. Mspm will be copy of spm.');
            fname = S.D;
            D = spm_eeg_load(fname);
            
            S2 = [];
            S2.montage.tra = eye(length(D.meegchannels));
            S2.montage.labelnew = D.chanlabels(D.meegchannels);
            S2.montage.labelorg = D.chanlabels(D.meegchannels);
            S2.D = D;
            S2.updatehistory = 0;
            S2.keepothers = 'yes';

            D = spm_eeg_montage(S2);
        end
        

    end
end
