function vbrt_infer_dynamic_model_gamma( X,options )

% vbrt_infer_dynamic_model_gamma( X,T,options )
    
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
    
    envset='tensorflow';
    
    disp(['Setting conda env: ' envset]);
    conda.setenv(envset)
    
    %%%%%%%
    % make the calls
    %python_cmd=['setup_tools.setup_dictionary(\"' options.workingdir '\", Q=' num2str(size(data,2)) ', ndicts=' num2str(size(data,2)) ', use_off_diags=False, use_greens_fns=False)'];  
    
    options.nsessions=1;
    options.vbrt.ntraining_init=100;
    options.vbrt.ntraining=options.vbrt.ntraining_init+200;
    options.vbrt.nportions=20;
    options.vbrt.subportion_length=140;
    
    options.vbrt.nstates=options.K;
    options.vbrt.model_mode='\"lstm\"';
    options.vbrt.model_name='\"dtfm\"';
    options.vbrt.beta_type='\"full\"';

    options.vbrt.do_recon='False';
    options.vbrt.do_fullprob_beta='True';
    options.vbrt.use_pca_cov_model='False';
    
    options.vbrt.npcs=40;
    
    options.vbrt.do_fullprob_theta='False';
    
    python_cmd=['infer_dynamic_model_gamma.reconstruct(\"' ...
        options.workingdir '\", \"' options.workingdir '\"' ...
        ', nsessions= ' num2str(options.nsessions) ...
        ', ntraining= ' num2str(options.vbrt.ntraining) ...
        ', ntraining_init= ' num2str(options.vbrt.ntraining_init) ...
        ', nportions=' num2str(options.vbrt.nportions) ...
        ', subportion_length=' num2str(options.vbrt.subportion_length) ...
        ', npcs=' num2str(options.vbrt.npcs) ...
        ', nstates=' num2str(options.vbrt.nstates) ...
        ', model_mode=' options.vbrt.model_mode ...
        ', load_model_epoch=None, epochs_per_model_save=100' ...
        ', model_name=' options.vbrt.model_name ...
        ', use_pca_cov_model=' options.vbrt.use_pca_cov_model ...
        ', beta_type=' options.vbrt.beta_type ...
        ', do_fullprob_theta=' options.vbrt.do_fullprob_theta ...
        ', do_fullprob_beta=' options.vbrt.do_fullprob_beta ...
        ', do_recon=' options.vbrt.do_recon ...
        ', state_model=' options.vbrt.state_model ...
        ')'];

    disp('Calling cmd:');
    cmd=['python /Users/woolrich/homedir/scripts/dynamic_network_recon_gamma/call_python_gamma.py --cmd="' python_cmd '"'];
    disp(cmd); 
    
    %runcmd(cmd);  
   
    %python_cmd=['plot_tools.plot_results(\"' ...
    %    options.workingdir '\", \"' options.workingdir '\"' ...
    %    ', time_range=[0,10000]' ...
    %    ')'];
       
    %disp('For python plots, call cmd:');
    %cmd=['python /Users/woolrich/homedir/scripts/dynamic_network_recon_gamma/call_python_gamma.py --cmd="' python_cmd '"'];
  
    %disp(cmd); 
    
    %runcmd(cmd);  
   
    %%

    %%%%%%%
    
    %pyversion(pyversion_old);
end

