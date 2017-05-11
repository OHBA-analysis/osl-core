function [bad_components, metrics, fig_handles] = identify_artefactual_components_manual(D,S)
    % AfRICA - ArteFact Rejection using Independent Component Analysis
    %
    % Syntax: [bad_components, fig_handles, fig_names, fig_titles] = identify_artefactual_components_manual(S)
    % S needs to contain be produced inside OSL_AFRICA.m
    %
    % OPTIONAL inputs to compute various metrics:
    %   -  ident.do_mains:       compute metric based on mains frequency [0 or 1]
    %   -  ident.do_kurt:        compute metric based on extreme kurtosis [0 or 1]
    %   -  ident.do_cardiac:     compute metric based on similarity to cardiac cycle [0 or 1]
    %   -  ident.artefact_channels: compute metric based on correlation with external channels
    %                               can be either part of a channame or a specific
    %                               channel idx [e.g. {'ECG','EOG','312'}]
    %   -  ident.launch_gui: if false, compute & save metrics & topos but don't use GUI to choose bad components
    % Output:
    %   - bad_components: A list of the components identified as bad.
    %
    % Romesh Abeysuriya 2017
    % Adam Baker 2014

    % Add defaults to ident params
    arg = inputParser;
    arg.KeepUnmatched = true; % Allow extra fields to be provided
    arg.addParameter('do_mains',false); 
    arg.addParameter('mains_frequency',50); 
    arg.addParameter('do_kurt',false); 
    arg.addParameter('do_cardiac',false); 
    arg.addParameter('launch_gui',true); 
    arg.addParameter('artefact_channels',{}); 
    arg.parse(S);
    S = arg.Results;

    % Compute metrics
    D = D.montage('switch',0);
    fig_handles = [];
    bad_components = []; 
    [metrics,tc] = compute_metrics(D,S);

    % Return just the components, or open the GUI
    if S.launch_gui
        bad_components = identify_artefactual_components_manual_gui(D,tc,D.ica.topos,metrics,D.ica.bad_components);
    end
