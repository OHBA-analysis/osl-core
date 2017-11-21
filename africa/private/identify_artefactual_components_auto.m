function [bad_components,reason] = identify_artefactual_components_auto(D,S)
    % Identify bad components automatically based on metrics computed from osl_africa() -> compute_metrics()
    % The settings struct S is the same one from osl_africa()
    %
    % Romesh Abeysuriya 2017
    % HL+MWW 2013

    metrics = D.ica.metrics;
    if any(D.ica.bad_components)
        fprintf(2,'Warning - bad components have already been marked.\nThese will be overwritten and replaced with the new run\n');
    end
    reason = cell(size(D.ica.sm,2),1); % Record reason for rejection

    % DETECT MAINS COMPONENTS
    if S.auto_do_mains
        mains_ind=find(metrics.mains.mains_frequency-1<metrics.mains.fmax(:) & metrics.mains.fmax(:) < metrics.mains.mains_frequency+1 & metrics.kurtosis.abs(:) < S.auto_mains_kurt_thresh);
        [~,idx] = sort(metrics.kurtosis.abs(mains_ind));
        mains_ind = mains_ind(idx);
        reason = reject(reason,mains_ind,S.auto_max_num_artefact_comps,'mains');
    end

    % Check correlations with any metrics whose name starts with
    % 'corr_chan'
    metric_names = fieldnames(metrics);
    corr_metrics = metric_names(cellfun(@(x) ~isempty(x),regexp(metric_names,'^corr_chan.*')));
        
    for j = 1:length(corr_metrics)
        to_reject = find(metrics.(corr_metrics{j}).value > S.auto_artefact_chans_corr_thresh);
        if isempty(to_reject)
            fprintf('No ICs correlated with Channel %s\n',strrep(corr_metrics{j},'corr_chan_',''));
        end
            
        for k = 1:length(to_reject)
            fprintf('Rejecting IC %d due to correlation with Channel %s (correlation = %.2f)\n',to_reject(k),strrep(corr_metrics{j},'corr_chan_',''),metrics.(corr_metrics{j}).value(to_reject(k)));
            reason{to_reject(k)} = sprintf('%s %s',reason{to_reject(k)},corr_metrics{j});
        end
    end
    
    % DETECT ARTEFACTS USING KURTOSIS   
    if S.auto_do_kurt
        [~,stats] = robustfit(ones(length(metrics.kurtosis.abs),1),metrics.kurtosis.abs,'bisquare',4.685,'off');
        outlier_inds = find(stats.w<S.auto_kurtosis_wthresh & metrics.kurtosis.abs>abs(S.auto_kurtosis_thresh) );
        [~,idx]=sort(stats.w(outlier_inds));
        outlier_inds = outlier_inds(idx);
        reason = reject(reason,outlier_inds,S.auto_max_num_artefact_comps,'kurtosis');
    end

    % ASSIGN THE BAD COMPONENTS
    bad_components = find(~cellfun(@isempty,reason));


function auto_reason = reject(auto_reason,candidates,limit,reason)
    % Take in a list of outliers identified for a reason, and a limit
    % Truncate and display message if required
    % Modify the artefact list and return it

    % Sort and truncate them
    if(length(candidates) > limit)
        fprintf(2,'%d components are kurtosis outliers. been detected. Will only label the worst %d components.',length(candidates),limit);
    end
    candidates = candidates(1:min(length(candidates),limit));

    if ~isempty(candidates)
        for j = 1:length(candidates)
            fprintf('Rejecting IC %d due to %s\n',candidates(j),reason);
            auto_reason{candidates(j)} = sprintf('%s %s',auto_reason{candidates(j)},reason);
        end
    else
        fprintf('No components rejected due to %s\n',reason);
    end

