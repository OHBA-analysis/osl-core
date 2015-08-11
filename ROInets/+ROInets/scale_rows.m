function y = scale_rows(x,s)
% SCALE_ROWS      Scale each row of a matrix.
% SCALE_ROWS(x,s) returns matrix y, same size as x, such that
% y(i,:) = s(i)*x(i,:)
% It is more efficient than diag(s)*x.

if isscalar(s),
    s = repmat(s, ROInets.rows(x), 1);
end%if

if iscolumn(s), 
    y = bsxfun(@times, x, s);
else
    y = bsxfun(@times, x, s');
end%if

% y = repmat(s(:), 1, ROInets.cols(x)).*x;
%y = (s(:)*ones(1,cols(x))).*x;
