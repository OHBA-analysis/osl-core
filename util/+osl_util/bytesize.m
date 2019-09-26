function [bsize,scale] = bytesize( b, unit )
% [bsize,scale] = osl_util.byte_size( b, unit='B' )
%
% Format input bytesize as a struct with conversions to kB, MB, GB, TB and PB
% If it is a BITsize (not BYTEsize), you can set the unit symbol to 'b' instead.
%
% If the input is not a scalar number, then we compute the bytesize of the argument 
% itself using whos.
%
% NOTE: Due to limitations of Matlab's whos function, bytesize estimation for 
% instances of handle classes is not accurate.
%
%
% The second output is the recommended unit to represent the input size.
% For example, the recommended unit for bytesize(123456) is 'kB'.
%
% If no output is collected, the size is printed to the console with the 
% recommended unit.
%
% JH

    if nargin < 2, unit='B'; end

    if ~isnumeric(b) || ~isscalar(b) 
        w = whos('b');
        b = w.bytes;
    end
    assert( isnumeric(b) && isscalar(b), 'Bad input bytesize.' );
    
    units = osl_util.mapfun( @(x) [x unit], {'','k','M','G','T','P'}, false );
    scale = units{min( numel(units), 1+floor( log(b)/log(1024) ) )};
    
    for i = 1:numel(units)
        u = units{i};
        bsize.(u) = b / 1024^(i-1);
    end
    
    if nargout == 0
        fprintf('%.2f %s\n', bsize.(scale), scale );
    end
    
end
