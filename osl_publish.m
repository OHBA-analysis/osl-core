function osl_publish(filename)
	% Make an HTML page on the OSL website
	%
	% INPUT - path of the file being published
	% OUTPUT - output html will be written in osl2/docs/matlab
	% 
	% To complete publishing, add files via Git and push to GitHub


	if nargin < 1 || isempty(filename) 
		filename = fullfile(getenv('OSLDIR'),'osl2','examples','publish_example.m');
	end
	
	publish(filename,'stylesheet',fullfile(getenv('OSLDIR'),'osl2','docs','mxdom2simplehtml_jekyll.xsl'),'format','html','outputDir',fullfile(getenv('OSLDIR'),'osl2','docs','matlab'))

