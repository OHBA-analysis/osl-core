function varargout = pfig(varargin)
    % Saves an image file using 'print'
    %
    % This is useful to quickly save things inside a loop. It provides more
    % control than saveas() because it allows the resolution in dpi to be set.
    % 300dpi is sufficient for a high resolution tiff for inclusion in papers
    %
    % Usage:
    %   pfig('test') - save current figure as 'test.png' to the desktop
    %   pfig('test.ext') - save current figure to the desktop with format specified by ext
    %   pfig(h,'test') - first argument can be a figure handle
    %   pfig(h,'test_%d',1) - additional arguments are used for string replacement

    % First, set the figure to have the same output dimension in the file
    % as it does on the screen
    
    % Romesh Abeysuriya 2017

    if nargin == 0
        if isempty(get(0,'Children')) % Exit immediately if no figures are open
            return
        end
        fhandle = gcf;
        fname = 'default.png';
    else
        
        if ishandle(varargin{1})
            fhandle = varargin{1};
            varargin = varargin(2:end);
        else
            fhandle = gcf;
        end

        fname = varargin{1};
     end

    % Perform optional string substitution
    if length(varargin) > 1
        fname = sprintf(fname,varargin{2:end});
    end

    [pathstr,fname,ext] = fileparts(fname);

    if isempty(pathstr)
        if isunix()
            pathstr = '~/Desktop';
        else
            pathstr = '%UserProfile%\Desktop';
        end
    end

    if isempty(ext)
        ext = '.png';
    end

    output_fname = fullfile(pathstr,[fname ext]);

    original_props = getprops(fhandle,{'PaperPositionMode','PaperUnits','Renderer','InvertHardcopy','PaperPosition','PaperSize','Color'});

    set(fhandle,'PaperPositionMode','auto');
    ppm = get(fhandle,'PaperPosition');
    set(fhandle,'PaperPosition',[0 0 ppm(3:4)]);
    set(fhandle,'PaperSize',ppm(3:4));

    if ~all(get(fhandle,'Color')==0.94) % The user set a custom background color
        set(fhandle,'InvertHardcopy','off'); % Use it
    elseif any(strcmp(ext,{'.eps','.pdf'}))
        set(fhandle,'InvertHardcopy','off'); % Use it
        set(fhandle,'Color','none');
    end

    switch ext
        case '.eps'
            print(fhandle,'-painters','-depsc',output_fname);
        case {'.jpg','.jpeg'} 
            print(fhandle,'-r300','-djpeg',output_fname);
        case {'.tif','.tiff'}
            print(fhandle,'-r300','-dtiff',output_fname);
        case {'.pdf'}
            print(fhandle,'-painters','-dpdf',output_fname);
        case {'.png'} 
            print(fhandle,'-r300','-dpng',output_fname);
        otherwise
            error(sprintf('Unknown extension format: %s',ext))
    end

    fprintf(1,'Image saved to %s\n',output_fname);

    set(fhandle,original_props);

    if nargout > 0
        varargout{1} = output_fname;
    end

    % This is convenient but system-dependent: removing background and trimming
    try
        if system('hash convert')==0 && strcmp(ext,'.png')
            system(sprintf('convert %s -trim -transparent white %s',output_fname,output_fname));
        end

        if system('hash pdfcrop')==0 && strcmp(ext,'.pdf')
            system(sprintf('pdfcrop %s ~/Desktop/temp.pdf >> /dev/null',output_fname));
            system(sprintf('mv ~/Desktop/temp.pdf %s',output_fname));
        end
    end


function s = getprops(h,proplist)
    s = struct;
    for j = 1:length(proplist)
        s.(proplist{j}) = get(h,proplist{j});
    end


