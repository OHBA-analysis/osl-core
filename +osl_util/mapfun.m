function out = mapfun( fun, val, unif )
%
% out = osl_util.mapfun( fun, val, unif=false )
%
% Use cellfun, arrayfun or structfun depending on the type of input.
% Returns a cell by default; set unif=true if you would like an array.
%
% It is fine if fun doesn't return anything, but then you should not collect the output.
%
% See also: cellfun, arrayfun, structfun
% 
% JH

    if nargin < 3, unif=false; end
    
    if iscell(val)
        map = @cellfun;
    elseif isscalar(val) && isstruct(val)
        map = @structfun;
    else
        map = @arrayfun;
    end
    
    if nargout == 0
        map( fun, val, 'UniformOutput', unif );
    else
        out = map( fun, val, 'UniformOutput', unif );
    end

end
