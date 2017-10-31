function artefacts(D)
	% Display bad stuff
	
% diagnostic plots
S=[];
S.D=D;
S.print_plots=1;
S.plotsdir=opt.results.plotsdir;
S.modality=opt.modalities{ii};
S.printprefix=[printprefix_mod '_outlier'];
S.plot_basename='double_maxfilter';
S.plot_basetitle='MAXFILTER: ';
[res fig_names fig_handles fig_titles]=osl_check_bad_chans(S);

report=osl_report_set_figs(report,fig_names,fig_handles,fig_titles);
report=osl_report_print_figs(report);





line 164 of detect_artefacts, plot stats.w etc
