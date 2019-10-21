function [ results ] = gumbel_infer_dynamic_model( X,options )

% [ results ] = gumbel_infer_dynamic_model( X,T,options )

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
    
    envset='tensorflow14';
    
    disp(['Setting conda env: ' envset]);
    conda.setenv(envset)
    
    %%%%%%%
    % make the calls
    %python_cmd=['setup_tools.setup_dictionary(\"' options.workingdir '\", Q=' num2str(size(data,2)) ', ndicts=' num2str(size(data,2)) ', use_off_diags=False, use_greens_fns=False)'];  
    
    options.vbrt=[];
    options.vbrt.ntraining_init=100;
    options.vbrt.ntraining=options.vbrt.ntraining_init+200;
    options.vbrt.nportions=20;
    options.vbrt.subportion_length=50;
    options.vbrt.nstates=options.K;
    options.vbrt.model_mode='\"lstm\"';
    options.vbrt.alpha_softxform_model='\"softmax\"';
    options.vbrt.model_name='\"dtfm\"';
    
    options.vbrt.do_recon='False';
    options.vbrt.do_fullprob_alpha='True';
    options.vbrt.do_fullprob_beta='True';
    options.vbrt.use_pca_cov_model='False';
    
    options.vbrt.npcs=20;
    
    options.vbrt.state_model='\"categorical\"';
    
    python_cmd=['infer_dynamic_model.reconstruct(\"' ...
        options.workingdir '\", \"' options.workingdir '\"' ...
        ', nsessions= ' num2str(length(fnames)) ...
        ', ntraining= ' num2str(options.vbrt.ntraining) ...
        ', ntraining_init= ' num2str(options.vbrt.ntraining_init) ...
        ', nportions=' num2str(options.vbrt.nportions) ...
        ', subportion_length=' num2str(options.vbrt.subportion_length) ...
        ', npcs=' num2str(options.vbrt.npcs) ...
        ', nstates=' num2str(options.vbrt.nstates) ...
        ', model_mode=' options.vbrt.model_mode ...
        ', alpha_softxform=' options.vbrt.alpha_softxform_model ...
        ', load_model_epoch=None, epochs_per_model_save=100' ...
        ', model_name=' options.vbrt.model_name ...
        ', use_pca_cov_model=' options.vbrt.use_pca_cov_model ...
        ', do_fullprob_alpha=' options.vbrt.do_fullprob_alpha ...
        ', do_fullprob_beta=' options.vbrt.do_fullprob_beta ...
        ', do_recon=' options.vbrt.do_recon ...
        ', state_model=' options.vbrt.state_model ...
        ')'];

    disp('Calling cmd:');
    disp(['python /Users/woolrich/homedir/scripts/dynamic_network_recon/call_python.py --cmd="' python_cmd '"']);
            
    runcmd(['python /Users/woolrich/homedir/scripts/dynamic_network_recon/call_python.py --cmd="' python_cmd '"'])   
    
    python_cmd=['plot_tools.plot_results(\"' ...
        options.workingdir '\", \"' options.workingdir '\"' ...
        ', time_range=[0,10000]' ...
        ')'];
       
    disp('Calling cmd:');
    disp(['python /Users/woolrich/homedir/scripts/dynamic_network_recon/call_python.py --cmd="' python_cmd '"']);
    
    runcmd(['python /Users/woolrich/homedir/scripts/dynamic_network_recon/call_python.py --cmd="' python_cmd '"'])   
    
    %%
    hmm=[];
    hmm.gamma=[];
    
    for ii=1:length(fnames)
        tmp=load([options.workingdir '/softxform_alpha_mean_store' num2str(ii-1) '.mat']);
        
        hmm.gamma=[hmm.gamma; tmp.softxform_alpha_mean_store0]
        
        keyboard;
        
        %hmm.statepath=[hmm.gamma; max(tmp.softxform_alpha_mean_store0)];
    end
        
    % options.K=8;workingdir  = [tilde '/homedir/vols_data/daisie/meg_data/'];i=1; 
    
    figure; 
    for kk=1:options.K, 
        subplot(options.K,1,kk);plot(hmm.gamma(:,kk),'LineWidth',1.5);
        a=axis;
        a(4)=1.1;
        axis(a);
    end;
    %%
    %%%%%%%
    
    %pyversion(pyversion_old);
end

