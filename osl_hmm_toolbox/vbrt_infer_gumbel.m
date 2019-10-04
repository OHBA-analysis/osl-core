function [ results ] = vbrt_infer_gumbel( X,options )

% [ results ] = vbrt_infer_dynamic_model_gamma( X,T,options )

    results = [];
    
    results.options=options;
    
    % save out each session in X to a directory
    clear fnames Hfnames;
    for ii=1:length(X)
        fnames{ii} = [options.workingdir, 'Y', num2str(ii)];
        data=X{ii}';
        data=normalise(data,2);
        save(fnames{ii}, 'data');
        
        Hfnames{ii} = [options.workingdir, 'H', num2str(ii)];
        nchans=size(data,1);
        
        H=eye(nchans);
        save(Hfnames{ii}, 'H');
    end
    
    % call python from matlab   
    
    %[~,pyversion_old]=pyversion();
    %pyversion('/Users/woolrich/anaconda/bin/python');
    
    envset='pytorch1';
    
    disp(['Setting conda env: ' envset]);
    conda.setenv(envset)
    
    cmd=['python /Users/woolrich/homedir/scripts/alex/test.py -d ' options.workingdir ' -n 1'];
    
    disp('Calling cmd:');
    disp(cmd); 
    
    %python_cmd=['plot_tools.plot_results(\"' ...
    %    options.workingdir '\", \"' options.workingdir '\"' ...
    %    ', time_range=[0,10000]' ...
    %    ')'];
       

    options.gumbel=[];
    
    options.data_path=options.workingdir;

    options.results_dirname=[options.workingdir 'results_test/'];
    

    cmd=['python /Users/woolrich/homedir/scripts/alex/model/src/gumbel.py ', ...
            '--data-dir ' options.data_path ' ', ...
            '--result-dir ' options.results_dirname ' ', ...
            '--temp 0.5 ' ...
            '--seq-length 400 ', ...
            '-n 800 ', ...
            '--batch-size 40 ', ...
            '--n-latent ' num2str(options.K) ' ', ...
            '--n-pc 0 ', ...
            '--n-hidden 32 ', ...
            '--n-layers 1 ', ...
            '--dropout 0.5 ', ...
            '--num-epochs 200 ', ...
            '--lr 0.001 ', ...
            '--rand-seed 42 ', ...
            '--train --infer ', ...
            '--checkpoint-every 20 ', ...
            '--bidirectional ', ...
            '--prior rnn ' ...
    ];
        
    disp('Calling cmd:');
    disp(cmd);
    
    %%
    
    results.gamma=[];
    
    for ii=1:length(fnames)
        
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
    tmp=load([options.workingdir  '/Djs_store.mat']);
    Djs=tmp.Djs_store;
    results.state=[];
    for kk=1:options.K
        results.state(kk).Omega.Gam_rate=squeeze(Djs(kk,:,:));
        ndim = length(results.state(kk).Omega.Gam_rate);
        results.state(kk).Omega.Gam_shape=ndim+2;
    end
    
    results.train=[];
    
    %%%%%%%
    
    %pyversion(pyversion_old);
end

