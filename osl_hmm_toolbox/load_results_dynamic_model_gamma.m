function results =load_results_dynamic_model_gamma(X,options)

    results = [];
    
    results.options=options;
    
    results.gamma=[];
    
    for ii=1:length(X)
        
        if strcmp(options.vbrt.state_model,'\"categorical\"')
            tmp=load([options.workingdir '/m_z_store' num2str(ii-1) '.mat']);
            alphas=tmp.m_z_store0;
        else
            tmp=load([options.workingdir '/m_alpha_store' num2str(ii-1) '.mat']);
            alphas=tmp.m_alpha_store0;
        end
        
        
        results.gamma=[results.gamma; alphas]
        
        %keyboard;
        
        %results.statepath=[results.gamma; max(tmp.softxform_alpha_mean_store0)];
    end
        
    % options.K=8;workingdir  = [tilde '/homedir/vols_data/daisie/meg_data/'];i=1; 
    
    %% plot state tcs
    if false
        figure; 
        for kk=1:options.K
            subplot(options.K,1,kk);plot(results.gamma(:,kk),'LineWidth',1.5);
            a=axis;
            a(4)=1.1;
            axis(a);
        end
    end
    
    %% load in and store state covariances
    tmp=load([options.workingdir  '/Djs_store.mat']);
    Djs=tmp.Djs_store;
    results.state=[];
    for kk=1:options.K
        results.state(kk).Omega.Gam_rate=squeeze(Djs(kk,:,:));
        ndim = length(results.state(kk).Omega.Gam_rate);
        results.state(kk).Omega.Gam_shape=ndim+2;
    end
    
    results.train=[];
    
    %%
    tmp=load([options.workingdir '/costs']);
    results.costs=tmp.costs;

end

