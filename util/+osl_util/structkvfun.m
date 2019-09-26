function out = structkvfun( fun, s, unif )
%
% out = osl_util.structkvfun( fun, s, unif=true )
%
% Variant of Matlab's structfun, which applies a function to each field of a struct-array.
% fun should be a function handle with TWO arguments:
%   fun( fieldname, fieldvalue )
%
% If unif=true (default), the output is an array with 
%   N rows = number of structures
%   F cols = number of fields. 
% For example, calling with a scalar struct yields a row-vector.
% Otherwise if unif=false, the output is a CELL.
%
% It is fine if fun does not return anything, but then you should not collect an output.
% 
% JH

    assert( isstruct(s), 'Second argument should be a structure.' );
    assert( isa(fun,'function_handle'), 'First argument should be a function handle.' );
    if nargin < 3, unif=true; end % array output by default

    n = numel(s);
    f = fieldnames(s);
    m = numel(f);

    if nargout > 0
        out = cell(n,m);
        for i = 1:n % structures
        for j = 1:m % fields
            out{i,j} = fun( f{j}, s(i).(f{j}) );
        end
        end
        if unif 
            out = cell2mat(out);
        end
    else
        for i = 1:n % structures
        for j = 1:m % fields
            fun( f{j}, s(i).(f{j}) );
        end
        end
    end

end
