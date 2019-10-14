function [ out ] = osl_tf_transform( S, dat )
%
% [ dattf ] = osl_tf_transform( S , dat )
%
% Either downsamples the time domain data
% or does a time-frequency transformation.
%
% If called with only S, and no dat, this function will only produce output
% containers and not do the transform.  If passed with a second argument
% dat, it will do the transformation on the data in dat.  dat should have
% dimensions [nTrials,nTimepoints]
%
% ARGUMENTS (in struct 'S')
% - S.tf_method                 : e.g. 'none', 'morlet', 'hilbert', 'hanning'
% - S.tf_freq_range             : the frequency window over which the tf is performed, e.g. [4 35]
% - S.tf_num_freqs                 : how many frequency bands to compute?
% - S.raw_times                  : the times in source_recon_results.times
% - S.tf_logtransform           : log transform the tf, trial by trial?
% - S.ds_factor                 : used for time domain, hilbert, morlet (a
%                                 value < 1 indicates lower temporal
%                                 resolution in TF output than was in the input)
% - S.tf_morlet_factor             : only used for morlet
% - S.tf_hanning_twin                   : used for hanning
% - S.tf_calc_amplitude
%
% IF S.ft_method = 'morlet', optional argument:
% - S.tf_morlet_basis (create using 'spm_eeg_morlet', for e.g.)
% >NB if this is NOT set, it will be reset each time the function is called
%
% OUTPUTS (in struct 'out')
% - out.tf_times
% - out.tf_freqs
% - out.tf_freq_res
%
% IF S.doTransform = 1 then additionally
% - out.dattf
%
% IF doing the HILBERT, then additionally
% - out.datbp
%
% IF doing the MORLET, then additionally
% - out.tf_morlet_basis

if nargin < 2
    S.doTransform = 0;
    dat = [];
else
    S.doTransform = 1;
end % if nargin < 2

if(~isfield(S,'tf_calc_amplitude')),
    S.tf_calc_amplitude=1;
end;

if(~S.tf_calc_amplitude && S.ds_factor~=1)
    %    error('Cannot downsample complex TF data');
end;

fsample=1/mode(diff(S.raw_times));

% Do downsampling / TF transform
switch S.tf_method

    case 'none'

        % always do
        out.tf_freqs = -1; % time domain
        if ~isempty(S.ds_factor) && S.ds_factor~=1,
            [dsn dsd]    = rat(S.ds_factor);
            tmp          = resample(S.raw_times,dsn,dsd);
            out.tf_times = linspace(S.raw_times(1),S.raw_times(end),numel(tmp));
            clear tmp
        else
            out.tf_times  = S.raw_times;
        end

        % Continue and do the transform if S.doTransform == 1
        if S.doTransform
            out.dattf = nan(size(dat,3),numel(out.tf_times),size(dat,1),numel(out.tf_freqs)); % [channel, time, trial, frequency]

            if ~isempty(S.ds_factor) && S.ds_factor~=1,
                [dsn dsd]    = rat(S.ds_factor);

                for iChan = 1:size(dat,1);

                    for iTrial = 1:size(dat,3); % indexes trials
                        dattrl=dat(iChan,:,iTrial);
                        dattrl=dattrl(:);
                        tmptrl    = resample(dattrl,dsn,dsd);
                        out.dattf(iChan,:,iTrial)=tmptrl;
                    end;
                end;

            else
                out.dattf = dat;
            end % if ~isempty(S.ds_factor)
        else

            % the first time the function is called, print some commentary
            if isempty(S.ds_factor) || S.ds_factor==1,
                tdHz  = fsample;
                glmHz=1/mode(diff(out.tf_times));
                disp(['TIME DOMAIN, NOT DOWNSAMPLING; raw data at ' num2str(tdHz) 'Hz, GLM at ' num2str(glmHz) 'Hz.'])
            else
                tdHz  = fsample;
                glmHz = 1/mode(diff(out.tf_times));
                disp(['TIME DOMAIN, DOWNSAMPLE FACTOR = ' num2str(S.ds_factor) '; raw data at ' num2str(tdHz) 'Hz, GLM at ' num2str(glmHz) 'Hz.']);
            end % if isempty(S.time_downsample_factor)

        end % if S.doTransform

    case 'hilbert'

        if ~isfield(S,'tf_hilbert_do_bandpass_for_single_freq')
            S.tf_hilbert_do_bandpass_for_single_freq=0;

        end;

        if(~S.tf_hilbert_do_bandpass_for_single_freq && S.tf_num_freqs == 1)
            if ~isempty(S.tf_freq_range)
                S.tf_hilbert_freq_ranges=S.tf_freq_range;
                % warning('S.tf_hilbert_do_bandpass_for_single_freq = 0. No further band-pass filtering is done during hilbert transform.  NB if your source recon band-pass is NOT THE SAME as the frequency band you want in the tf-transform, you need to set oat.first_level.tf_hilbert_do_bandpass_for_single_freq = 1');
            end;
        end;

        % always do
        if(~isempty(S.tf_hilbert_freq_ranges)),
            if S.tf_num_freqs~=size(S.tf_hilbert_freq_ranges,1)
                error('S.tf_num_freqs~=size(S.tf_hilbert_freq_ranges,1)');
            end;
            out.tf_freq_res = [];
            out.tf_freq_ranges = S.tf_hilbert_freq_ranges;
            out.tf_freqs   = mean(out.tf_freq_ranges,2)';
        elseif S.tf_num_freqs > 1,
            out.tf_freq_res = S.tf_hilbert_freq_res;
            out.tf_freqs = linspace(S.tf_freq_range(1)+out.tf_freq_res/2, S.tf_freq_range(2)-out.tf_freq_res/2, S.tf_num_freqs);

            for f = 1 : S.tf_num_freqs
                out.tf_freq_ranges(f,:)=[out.tf_freqs(f)-out.tf_freq_res/2, out.tf_freqs(f)+out.tf_freq_res/2];
            end;
        else
            out.tf_freq_res = S.tf_hilbert_freq_res;
            out.tf_freqs   = (S.tf_freq_range(2) + S.tf_freq_range(1)) ./ 2;

            for f = 1 : S.tf_num_freqs
                out.tf_freq_ranges(f,:)=[out.tf_freqs(f)-out.tf_freq_res/2, out.tf_freqs(f)+out.tf_freq_res/2];
            end;
        end % if S.tf_num_freqs > 1

        out.tf_freq_ranges=max(1e-3,out.tf_freq_ranges);
        out.tf_freqs   = mean(out.tf_freq_ranges,2)';

        % downsample the time vector appropriately
        if ~isempty(S.ds_factor) && S.ds_factor~=1,
            [dsn dsd]    = rat(S.ds_factor);
            tmp          = resample(S.raw_times,dsn,dsd);
            out.tf_times = linspace(S.raw_times(1),S.raw_times(end),numel(tmp));
            clear tmp
        else
            out.tf_times  = S.raw_times;
        end

        dat = permute(dat,[3,2,1]);

        freq_ind=cell(length(out.tf_freqs),1);
        % Continue and do the transform if S.doTransform == 1
        if S.doTransform
            ft_progress('init', 'text', 'Doing Hilbert transform...')      % ascii progress bar
            out.dattf = nan(size(dat,3),numel(out.tf_times),size(dat,1),numel(out.tf_freqs)); % [channel, time, trial, frequency]
            out.datbp = nan(size(dat,3),numel(out.tf_times),size(dat,1),numel(out.tf_freqs)); % [channel, time, trial, frequency]
            for iChan = 1:size(dat,3);

                ft_progress(iChan/size(dat,3));
                tempoutpow = zeros(size(dat,1),length(S.raw_times),length(out.tf_freqs));
                tempoutbp = zeros(size(dat,1),length(S.raw_times),length(out.tf_freqs));

                for iTrial = 1:size(dat,1) % indexes trials
                    dattrl=dat(iTrial,:,iChan);
                    dattrl=dattrl(:);

                    for f = 1 : length(out.tf_freqs)
                        if(~S.tf_hilbert_do_bandpass_for_single_freq && S.tf_num_freqs == 1)
                            tempdatbp=dattrl;

                        else
                            % bandpass data
                            bp_cutoff=out.tf_freq_ranges(f,:);
                            [tempdatbp freq_ind{f}]=bandpass(dattrl,bp_cutoff,fsample,0,freq_ind{f});
                        end;

                        signal_h = hilbert(tempdatbp);

                        % square the hilbert transform and take the
                        % +ve square root to get the estimated envelope
                        if(S.tf_calc_amplitude)
                            tempdatpow = (sqrt(signal_h.*conj(signal_h)));
                        else
                            tempdatpow=signal_h;
                        end;

                        tempoutpow(iTrial,:,f) = tempdatpow;
                        tempoutbp(iTrial,:,f)  = tempdatbp;

                    end % f = 1 : length(out.tf_freqs)
                end % for iTrial = 1:size(dat,1);

                % if requested, downsample
                if ~isempty(S.ds_factor) && S.ds_factor~=1,
                    [dsn dsd]    = rat(S.ds_factor);

                    nfreqs = size(tempoutpow,3);

                    tempoutpow   = permute(tempoutpow,[2,1,3]);
                    tempoutbp    = permute(tempoutbp,[2,1,3]);

                    tempoutpow   = reshape(tempoutpow,[size(tempoutpow,1),(size(tempoutpow,2)*nfreqs)]);
                    tempoutbp    = reshape(tempoutbp, [size(tempoutbp,1), (size(tempoutbp, 2)*nfreqs)]);

                    tempoutpow   = resample(tempoutpow,dsn,dsd);
                    tempoutbp    = resample(tempoutbp,dsn,dsd);

                    tempoutpow   = reshape(tempoutpow,[size(tempoutpow,1),(size(tempoutpow,2)./nfreqs),nfreqs]);
                    tempoutbp    = reshape(tempoutbp,[size(tempoutbp,1), (size(tempoutbp,2)./nfreqs),nfreqs]);

                    out.dattf(iChan,:,:,:) = permute(tempoutpow,[4,1,2,3]); % [channel, time, trial, frequency]
                    clear tempoutpow
                    out.datbp(iChan,:,:,:) = permute(tempoutbp,[4,1,2,3]);
                    clear tempoutbp

                else
                    out.dattf(iChan,:,:,:) = permute(tempoutpow,[4,2,1,3]); % [channel, time, trial, frequency]
                    clear tempoutpow
                    out.datbp(iChan,:,:,:) = permute(tempoutbp,[4,2,1,3]);
                    clear tempoutbp
                end
            end % for iChan
            ft_progress('close')
        else

            % the first time the function is called, print some commentary
            disp('USING HILBERT TF');

            disp(['Freq ranges: ' mat2str(out.tf_freq_ranges,3)]);
            disp(['Freq resolution: ' num2str(out.tf_freq_res,3)]);

        end % if S.doTransform


    case 'morlet'

        % always do
        if S.tf_num_freqs == 1
            out.tf_freqs   = (S.tf_freq_range(1) + S.tf_freq_range(2))/2;
            out.tf_freq_res =  S.tf_freq_range(1) - S.tf_freq_range(2);
        else
            out.tf_freqs   = linspace(S.tf_freq_range(1),S.tf_freq_range(2),S.tf_num_freqs);
            out.tf_freq_res = out.tf_freqs(2) - out.tf_freqs(1);
        end % if S.tf_num_freqs == 1

        % downsample the time vector appropriately
        if ~isempty(S.ds_factor) && S.ds_factor~=1,
            [dsn dsd]    = rat(S.ds_factor);
            tmp          = resample(S.raw_times,dsn,dsd);
            out.tf_times = linspace(S.raw_times(1),S.raw_times(end),numel(tmp));
            clear tmp
        else
            out.tf_times  = S.raw_times;
        end % if ~isempty(S.ds_factor)

        if ~isfield(S,'tf_morlet_basis')
            disp('Creating Morlet basis set.  If you are seeing message many times, you may wish to pass a morlet basis set to ''osl_tf_transform''');
            fres=fsample;
            out.tf_morlet_basis = spm_eeg_morlet(S.tf_morlet_factor, 1000/fres, out.tf_freqs); 
            S.tf_morlet_basis   = out.tf_morlet_basis;
        end % if ~isfield(S,'tf_morlet_basis')

        % Continue and do the transform if S.doTransform == 1
        if S.doTransform

            dat = permute(dat,[3,2,1]);

            ft_progress('init', 'text', 'Doing Morlet transform...')      % ascii progress bar
            out.dattf = nan(size(dat,3),numel(out.tf_times),size(dat,1),numel(out.tf_freqs)); % [channel, time, trial, frequency]
            out.datbp = nan(size(dat,3),numel(out.tf_times),size(dat,1),numel(out.tf_freqs)); % [channel, time, trial, frequency]

            for iChan = 1:size(dat,3);

                ft_progress(iChan/size(dat,3));

                tempdattf = zeros(size(dat,1),length(S.raw_times), length(out.tf_freqs));
                for iTrial = 1:size(dat,1); % indexes trials

                    dattrl = dat(iTrial,:,iChan);
                    dattrl = dattrl(:);

                    for f = 1 : length(out.tf_freqs)
                        tmp = conv(dattrl, S.tf_morlet_basis{f},'same');

                        % power
                        if(S.tf_calc_amplitude)
                            tmp = sqrt(tmp.*conj(tmp))'; % tf_num_freqs x num_timepoint
                        end;

                        tempdattf(iTrial,:,f) = tmp;
                    end % for f = 1 : length(out.tf_freqs)
                end % for iTrial = 1:size(dat,1);

                % downsample
                if ~isempty(S.ds_factor) && S.ds_factor~=1,
                    [dsn dsd]    = rat(S.ds_factor);

                    nfreqs = size(tempdattf,3);

                    tempdattf   = permute(tempdattf,[2,1,3]);
                    tempdattf   = reshape(tempdattf,[size(tempdattf,1),(size(tempdattf,2)*nfreqs)]);
                    tempdattf   = resample(tempdattf,dsn,dsd);
                    tempdattf   = reshape(tempdattf,[size(tempdattf,1),(size(tempdattf,2)./nfreqs),nfreqs]);

                    out.dattf(iChan,:,:,:) = permute(tempdattf,[4,1,2,3]); % [channel, time, trial, frequency]
                    clear tempdattf;

                else
                    tempdattf   = permute(tempdattf,[2,1,3]);
                    out.dattf(iChan,:,:,:) = permute(tempdattf,[4,1,2,3]); % [channel, time, trial, frequency]
                    clear tempdattf;
                end
            end % for iChan
            ft_progress('close')
        else

            % the first time the function is called, print some commentary
            disp('USING MORLET TF');

            disp(['Freqs: ' num2str(out.tf_freqs)]);

        end % if S.doTransform

    case 'hanning'

        fs = 1./mode(diff(S.raw_times));
        tf_hanning_timestep = 1 ./ (S.ds_factor * fs);

        % We need to set the frequency bands such that we have a full number
        % of cycles within the time window.
        if S.tf_num_freqs > 1
            out.tf_freqs   = linspace(S.tf_freq_range(1),S.tf_freq_range(2),S.tf_num_freqs);
        else
            out.tf_freqs = mean([S.tf_freq_range(1),S.tf_freq_range(2)]);
        end
        cycles = S.tf_hanning_ncycles;


        S.tf_hanning_biggest_twin = cycles.*(1./out.tf_freqs(1));

        if S.tf_hanning_biggest_twin > (S.raw_times(end) - S.raw_times(1))
            error('The maximum time window specified for the analysis is less than one cycle of the lowest frequency!')
        end

        % determine the true time range we can use
        first_viable_time = S.raw_times(1)   + S.tf_hanning_biggest_twin/2 + tf_hanning_timestep;
        last_viable_time  = S.raw_times(end) - S.tf_hanning_biggest_twin/2 - tf_hanning_timestep;
        out.tf_times    = first_viable_time:tf_hanning_timestep:last_viable_time;

        if S.doTransform % ... continue and do the transform

            % parse data matrix into fieldtrip format: cell array
            % of trials; as this is voxelwise, each array will be
            % {1xn}
            ftraw = [];
            ftraw.trial = cell(size(dat,3),1);
            for iTrial = 1:size(dat,3);
                ftraw.trial{iTrial} = dat(:,:,iTrial);
            end % for iTrial = 1:size(dat,1);
            ftraw.time = {S.raw_times};
            ftraw.time = repmat(ftraw.time,[size(dat,3),1]);
            ftraw.fsample = fsample;
            for ichan = 1:size(dat,1)
                ftraw.label{ichan, 1} = ['Ch' num2str(ichan)];
            end

            cfg = [];
            cfg.output      ='fourier'; % NB could be 'fourier' - may need care though...
            cfg.taper       ='hanning';
            cfg.method      ='mtmconvol';
            cfg.foi         = out.tf_freqs;
            cfg.toi         = out.tf_times;
            cfg.t_ftimwin   = cycles./cfg.foi;
            cfg.keeptrials = 'yes';

            freq = ft_freqanalysis(cfg, ftraw);
            out.dattf = permute(freq.fourierspctrm,[2 4 1 3]); % [channel,time,trial,frequency]
        else

            disp(['Hanning taper TF.']);
            disp(['The first time for which the TF will be calculated will be ' num2str(first_viable_time) 's.']);
            disp(['The last time for which the TF will be calculated will be ' num2str(out.tf_times(end)) 's.']);

        end % if S.doTransform
end % switch S.tf_method
