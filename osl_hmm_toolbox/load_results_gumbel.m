function results=load_results_gumbel(X,options)

    results = [];
    
    results.options=options;
    
    results.gamma=[];
    
    for ii=1:length(X)
        
        tmp=readNPY([options.workingdir 'results_test/q_z_all.npy']);
        
        alphas=reshape(tmp,[size(tmp,1)*size(tmp,2),size(tmp,3)]);
        
        results.gamma=[results.gamma; alphas];
        
        %keyboard;
        
    end
    
    %%
    tmp=readNPY([options.workingdir 'results_test/train_loss.npy']);
    results.costs=tmp;

    
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
    tmp=readNPY([options.workingdir 'results_test/inferred_covariances.npy']);
    Djs=tmp;
    results.state=[];
    for kk=1:options.K
        results.state(kk).Omega.Gam_rate=squeeze(Djs(kk,:,:));
        ndim = length(results.state(kk).Omega.Gam_rate);
        results.state(kk).Omega.Gam_shape=ndim+2;
    end
    
    results.train=[];

