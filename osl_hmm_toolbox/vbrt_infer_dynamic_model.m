function [ results ] = vbrt_infer_dynamic_model( X,options )

% [ results ] = vbrt_infer_dynamic_model( X,T,options )

    results = [];
    
    results.options=options;
    
    % save out each session in X to a directory
    for ii=1:length(X)
        fname = [options.workingdir, 'X', num2str(ii)];
        data=X{ii};
        save(fname, 'data');
    end
    
    % call python from matlab   
    
    %[~,pyversion_old]=pyversion();
    %pyversion('/Users/woolrich/anaconda/bin/python');
    
    conda.setenv('tensorflow')
    
    %%%%%%%
    % make the calls
    python_cmd=['setup_tools.setup_dictionary(\"' options.workingdir '\", Q=' num2str(size(data,2)) ', ndicts=' num2str(size(data,2)) ', use_off_diags=False, use_greens_fns=False)'];  
    runcmd(['python /Users/woolrich/Dropbox/vols_scripts/dynamic_network_recon/dynamic_network_recon/call_python.py --cmd="' python_cmd '"'])
    
    
    %%%%%%%
    
    %pyversion(pyversion_old);
end

