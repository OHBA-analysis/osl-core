function osl_publish(filename)
	% Provide a filename to publish a specific file
	% or leave empty to regenerate everything
	if nargin < 1 || isempty(filename) 
		filename = fullfile(getenv('OSLDIR'),'osl2','examples','publish_example.m');
	end
	
	%publish(filename,'format','html','outputDir',fullfile(getenv('OSLDIR'),'osl2','docs','_matlab'))
	publish(filename,'stylesheet',fullfile(getenv('OSLDIR'),'osl2','docs','mxdom2simplehtml_jekyll.xsl'),'format','html','outputDir',fullfile(getenv('OSLDIR'),'osl2','docs','matlab'))

