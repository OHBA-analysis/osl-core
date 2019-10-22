function s = repstruct( fields, varargin )
%
% s = osl_util.repstruct( fields, varargin )
% s = osl_util.repstruct( fields, varargin )
%
% Create a struct-array with specified fields.
% This is equivalent to repmat( struct('field1',[],'field2',[],...), varargin{:} ).
% 
% This function can also be called to create an "empty" struct, where fieldnames are specified, 
% but values are not. For example:
%   s = osl_util.repstruct( {'a','b','c'} )
% creates a scalar structure with fields {'a','b','c'}.
%
% See also: repmat
%
% JH

    assert( iscellstr(fields), 'Expects a cell of fieldnames in input.' );
    if nargin < 2
        args = {1};
    else
        args = varargin;
    end

    n = numel(fields);
    f = cell(1,2*n);
    f(1:2:end) = fields;

    s = repmat( struct(f{:}), args{:} );

end
