function varargout = runcmd(cmd,varargin)
    % Wrapper to capture output of system calls. Supports sprintf token replacement
    % 
    % [output,return_value] = runcmd(cmd,varargin)
    %
    % INPUTS
    % - cmd - command to run
    % - varargin - String substitutions for cmd (same syntax as sprintf)
    % 
    % OUTPUTS
    % output - terminal output from command
    % return_value - value returned by the command
    % 
    % EXAMPLE USAGE
    % - runcmd('ls') % Run a command, print output to terminal
    % - output = runcmd('ls') % Run command and capture output strings
    % - output = runcmd('ls %s',getenv('OSLDIR')) % Use string substitution


    if nargin > 1
        cmd = sprintf(cmd,varargin{:});
    end

    [return_value, output] = system(cmd);

    if(return_value ~= 0)
        throw(MException('runcmd:error',sprintf('runcmd call:\n%s\nReturn value: %d\nProduced error:\n%s\n',cmd,return_value,output)));
    end

    % If no output is captured and the function produced no output, don't print anything to the terminal
    if nargout > 0 || ~isempty(output)
        varargout{1} = output;
    end

    if nargout > 1
        varargout{2} = return_value;
    end