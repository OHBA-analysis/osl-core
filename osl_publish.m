function osl_publish(filename,evalCode)
	% Make an HTML page on the OSL website
	%
	% INPUT - path of the file being published
	% OUTPUT - output html will be written in osl2/docs/matlab
	% 
	% To complete publishing, add files via Git and push to GitHub

	% TODO - see if it ends up being just stuff in a couple of folders that gets published - if so, automatically build all of them
	% TODO - use 'git hash-object' to see if a file needs to be republished
	% TODO - maybe do it automatically?
	
	if nargin < 2 || isempty(evalCode) 
		evalCode = true;
	end
	
	if nargin < 1 || isempty(filename) 
		filename = fullfile(osldir,'osl2','examples','publish_example.m');
	end
	
	output_html = publish(filename,'evalCode',evalCode,'stylesheet',fullfile(osldir,'osl2','docs','mxdom2simplehtml_jekyll.xsl'),'format','html','outputDir',fullfile(osldir,'osl2','docs','matlab'));

	runcmd('open %s',output_html)