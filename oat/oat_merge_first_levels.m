
function [oat] = firstLevelSessCombine(sessionInds,oat)
%
% [oat] = firstLevelSessCombine(sessionInds,oat)
%
% Takes a cell array of first-level result file indices grouped by subject,
% i.e.
% {[1,2],[3,4],...etc} if 2 sessions per subject, and combines the
% session-wise first level GLMs into a subject-wise first level GLM using
% flame1

% note that resultscope here are:
% num_tpts x  num_firstlevel_contrasts x num_freqs x num_sessions
% ...and the GLM fitted by flame1 is computing an average over sessions

results_fnames = cell(numel(sessionInds),1);
for iSub = 1:numel(sessionInds)
   
    xnew = [];
    % get the first_level data
    for iSess = 1:numel(sessionInds{iSub})
        fname = oat.first_level.results_fnames{sessionInds{iSub}(iSess)};
        tmp   = load([oat.source_recon.dirname '/' fname]);
        res_strct(iSess).res   = tmp.oat_stage_results;
       
        xnew = [xnew ; res_strct(iSess).res.x];
       
        % delete the session-wise first level files once the data is
        % extracted from them
        delete([oat.source_recon.dirname '/' fname '.mat']);
       
    end % for iSess = 1:numel(sessionIndices{iSub})
   
   
    % get the number of voxels, times, freqs, sessions
    nVox        = numel(res_strct(1).res.mask_indices_in_source_recon);
    nTimes      = numel(res_strct(1).res.times);
    nFreqs      = numel(res_strct(1).res.frequencies);
    nContrasts  = numel(res_strct(1).res.contrasts);
    nSessions   = numel(sessionInds{iSub});
   
    % make containers
    first_level_results.cope    = nan(nVox,nTimes,nContrasts,nFreqs);
    first_level_results.stdcope = nan(nVox,nTimes,nContrasts,nFreqs);
   
    % voxel loop
    ft_progress('init', 'etf');
    for iVox = 1:nVox
        ft_progress(iVox/nVox);
        % time loop
        for iTime = 1:nTimes
            % contrast loop
            for iContrast = 1:nContrasts
                % freq loop
                for iFreq = 1:nFreqs
                   
                    % session loop to concatenate sessions
                    resultscope   = nan(nTimes,nContrasts,nFreqs,nSessions);
                    resultsstdcope = nan(nTimes,nContrasts,nFreqs,nSessions);
                    for iSess = 1:numel(sessionInds{iSub})
                        if nFreqs == 1;
                            resultscope(iTime,iContrast,iFreq,iSess)    = res_strct(iSess).res.cope(iVox,iTime,iContrast);
                            resultsstdcope(iTime,iContrast,iFreq,iSess) = res_strct(iSess).res.stdcope(iVox,iTime,iContrast);
                        else
                            resultscope(iTime,iContrast,iFreq,iSess)    = res_strct(iSess).res.cope(iVox,iTime,iContrast,iFreq);
                            resultsstdcope(iTime,iContrast,iFreq,iSess) = res_strct(iSess).res.stdcope(iVox,iTime,iContrast,iFreq);
                        end
                    end % for iSess = 1:numel(sessionInds{iSub})
                   
                    cope=permute(resultscope(iTime,iContrast,iFreq,:),[4 1 2 3]);
                    stdcope=permute(resultsstdcope(iTime,iContrast,iFreq,:),[4 1 2 3]);
                   
                    S=stdcope.^2;
                    z=ones(size(S,1),1);
                    [gam, ~, covgam] = flame1(cope,z,S,1);
                   
                    first_level_results.cope(iVox,iTime,iContrast,iFreq)=gam;
                    first_level_results.stdcope(iVox,iTime,iContrast,iFreq)=sqrt(covgam);
                   
                end % for iFreq = 1:nFreqs
            end % for iContrast = 1:nContrasts
        end % for iTime = 1:nTimes
    end % for iVox = 1:nVox
    ft_progress('close');
   
    % save the new first level results
    % first_level_results.cope already specified
    % first_level_results.stdcope already specified
    first_level_results.mask_indices_in_source_recon = res_strct(1).res.mask_indices_in_source_recon; % does this differ between sensor and source space, or can it be set the same way here?
    first_level_results.times               = res_strct(1).res.times;
    first_level_results.frequencies         = res_strct(1).res.frequencies;
    % (stdcope,cope, set above)
    first_level_results.pseudo_zstat_var    = res_strct(1).res.psudo_zstat_var; % probably wrong?
    first_level_results.D_sensor_data       = res_strct(1).res.D_sensor_data; % almost certainly wrong?
    first_level_results.chanind             = res_strct(1).res.chanind; % probably only for a sensor level analysis - wrap in a conditional?
    first_level_results.x                   = xnew;
    first_level_results.contrasts           = res_strct(1).res.contrasts;
    first_level_results.source_recon_results_fname = res_strct(1).res.source_recon_results_fname;
    first_level_results.S                   = res_strct(1).res.S;
    first_level_results.bf_S                = res_strct(1).res.bf_S;
    first_level_results.level               = res_strct(1).res.level;
    first_level_results.name                = res_strct(1).res.name;
    first_level_results.recon_method        = res_strct(1).res.recon_method;
   
    % will gridstep need to be specified for source space?
    if isfield(res_strct(1).res,'mni_coord') % is this still a field for beamformed data?
        first_level_results.mni_coord = res_strct(1).res.mni_coord;
    end
   
    first_level_results.subject_name = ['subject' num2str(iSub)];
    first_level_results.fname=[ first_level_results.subject_name '_' first_level_results.name ];
   
    disp(['Saving session-concatenated statistics in file ' oat.name '/' first_level_results.fname]);
    oat_save_results( oat, first_level_results );
    results_fnames{iSub}=first_level_results.fname;
   
end % for iSub = 1:numel(sessionInds)   

% update oat.first_level
oat.first_level.results_fnames = results_fnames;
