function oat_save_results( oat, oat_stage_results )

%oat_save_results( oat, oat_stage_results )

oat_stage_results.osl2_version=osl2_version;

save([oat.source_recon.dirname '/' oat_stage_results.fname], 'oat_stage_results', '-v7.3');
