function osl_publish(filename,evalCode,stylesheet)
	% Make an HTML page on the OSL website
	%
	% INPUTS
	% - filename - name of the file being published
	% - evalCode - run the code and save code output and figures (default=true)
	% - stylesheet - specify which formatting to apply (default=github pages formatting, set to '' for default matlab formatting)
	%
	% OUTPUT - output html will be written in osl2/docs/matlab
	% 
	% To complete publishing, add files via Git and push to GitHub
	%
	% PUBLISHING WORKFLOW
	% 1. After writing your script, use osl_publish('myscript',false,'') to do an initial test of the formatting
	% 2. Then, use osl_publish('myscript',true,'') to check the figures are correct
	% 3. Finally, use osl_publish('myscript') to generate the github-pages ready files
	% 4. Add the new files in the docs folder, commit, and push, to publish them online
	%
	%
	% Romesh Abeysuriya 2017
	
	if nargin < 3
		stylesheet = fullfile(osldir,'osl2','docs','mxdom2simplehtml_jekyll.xsl');
	end
	
	if nargin < 2 || isempty(evalCode) 
		evalCode = true;
	end
		
	output_html = publish(filename,'evalCode',evalCode,'stylesheet',stylesheet,'format','html','outputDir',fullfile(osldir,'osl2','docs','matlab'));

	runcmd('open %s',output_html)