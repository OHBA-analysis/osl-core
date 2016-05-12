function opt_save_results( opt, opt_results )

% opt_save_results( opt, opt_results )
%
% saves opt results for a single session

opt_results.osl2_version=osl2_version;
save([opt.dirname '/' opt_results.fname], 'opt_results', '-v7.3');
