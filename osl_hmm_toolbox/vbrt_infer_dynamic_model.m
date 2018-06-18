function [ results ] = vbrt_infer_dynamic_model( X,options )

% [ results ] = vbrt_infer_dynamic_model( X,T,options )

    results = [];
    
    results.options=options;
    
    % save out each session in X to a directory
    clear fnames Hfnames;
    for ii=1:length(X)
        fnames{ii} = [options.workingdir, 'Y', num2str(ii)];
        data=X{ii}';
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
    conda.setenv('tensorflow14')
    
    %%%%%%%
    % make the calls
    %python_cmd=['setup_tools.setup_dictionary(\"' options.workingdir '\", Q=' num2str(size(data,2)) ', ndicts=' num2str(size(data,2)) ', use_off_diags=False, use_greens_fns=False)'];  
    
    options.vbrt=[];
    options.vbrt.n_training=400;
    options.vbrt.n_portions=10;
    options.vbrt.subportion_length=400;
    options.vbrt.npcs=10;
    options.vbrt.nfactors=options.K;
    options.vbrt.model_mode='\"lstm\"';
    options.vbrt.alpha_softxform_model='\"softmax\"';
    options.vbrt.model_name='\"dtfm\"';
    
    python_cmd=['infer_dynamic_model.reconstruct(\"' ...
        options.workingdir '\", \"' options.workingdir '\"' ...
        ', n_sessions= ' num2str(length(fnames)) ...
        ', n_training= ' num2str(options.vbrt.n_training) ...
        ', n_portions=' num2str(options.vbrt.n_portions) ...
        ', subportion_length=' num2str(options.vbrt.subportion_length) ...
        ', npcs=' num2str(options.vbrt.npcs) ...
        ', nfactors=' num2str(options.vbrt.nfactors) ...
        ', model_mode=' options.vbrt.model_mode ...
        ', alpha_softxform=' options.vbrt.alpha_softxform_model ...
        ', load_model_epoch=None, epochs_per_model_save=100' ...
        ', model_name=' options.vbrt.model_name ')' ];
     
    disp('Calling cmd:');
    disp(['python /Users/woolrich/homedir/scripts/dynamic_network_recon/call_python.py --cmd="' python_cmd '"']);
    
    runcmd(['python /Users/woolrich/homedir/scripts/dynamic_network_recon/call_python.py --cmd="' python_cmd '"'])
    
    
    %%%%%%%
    
    %pyversion(pyversion_old);
end

