function output_html = osl_publish( filename, evalCode, matlabStyle, docsFolder )
%
% Make an HTML page on the OSL website
%
% INPUTS
% - filename - name of the file being published
% - evalCode - run the code and save code output and figures (default=true)
% - docsFolder - specify the root of the osl-docs repo (default=osldir('osl-docs'))
% - matlabStyle - set to true to disable GitHub pages formatting (default=false)
%
% OUTPUT - output html will be written in osl-docs/matlab
%
% To complete publishing, add files via Git and push to GitHub
%
%
% PUBLISHING WORKFLOW
%
% 1. After writing your script, use osl_publish('myscript',false) to do an initial test of the formatting
% 2. Then, use osl_publish('myscript',true) to check the figures are correct
% 3. Finally, use osl_publish('myscript',true) to generate the github-pages ready files
% 4. Add the new files in the docs folder, commit, and push, to publish them online
%
% NOTE: if the script to be published is inside a folder added to the path (e.g. example/osl_example_blah),
%       you do NOT need to specify that folder when calling osl_publish. That is, calling
%           osl_publish( 'osl_example_blah', true )
%       will work just fine.
%
% RA 2017, JH 2019

    if nargin < 4 || isempty(docsFolder)
        docsFolder = osldir('osl-docs');
    end

    if nargin < 3 || isempty(matlabStyle) 
        matlabStyle = false; 
    end
    
    if nargin < 2 || isempty(evalCode)
        evalCode = true;
    end
    
    % set options
    opt.format = 'html';
    opt.evalCode = evalCode;
    opt.outputDir = fullfile(docsFolder,'matlab');
    
    if ~matlabStyle
        opt.stylesheet = fullfile( docsFolder, 'mxdom2simplehtml_jekyll.xsl' );
    end

    % publish
    output_html = publish( filename, opt );

    % Command only works on OSX. (Linux could use xdg-open instead.)
    % runcmd('open %s',output_html)
    
end