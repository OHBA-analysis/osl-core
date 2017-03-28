function [results_fnames results] = osl_run_first_level( oat )

% [results_fnames results] = osl_run_first_level( oat )
%
% takes in an OAT, which needs to be setup by calling oat=osl_setup_oat(S), struct
% and runs first level analysis
%
% This function should normally be called using osl_run_oat(oat);
%
% MWW 2011

OSLDIR = getenv('OSLDIR');

if ~oat.first_level.is_epoched    

    switch lower(oat.source_recon.type)
        case 'scalar'
           
                [results_fnames results] = oat_run_first_level_continuous(oat); % Do time-wise analysis state-wise        

        otherwise   
                
            error('Not supported');
    end

else
    
    [results_fnames results] = oat_run_first_level_epoched(oat); % Do time-wise analysis state-wise       
    
end
