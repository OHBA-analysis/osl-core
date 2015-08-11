function z = p_to_z_two_tailed(p, sign_z)
%Z_TO_P_TWO_TAILED convert p-value to standard z-value in two-tailed test
%
% Z = Z_TO_P_TWO_TAILED(P, SIGN_Z) converts P to standard Z-values,
%   interpreting P as a two-tailed p-value. Z values are returned with
%   signs provided in SIGN_Z. 

if nargin < 2 || isempty(sign_z), 
    sign_z = 1; 
end%if

assert(isscalar(sign_z) || all(size(p) == size(sign_z)), ...
       [mfilename ':IncompatibleInputSizes'], ...
       'SIGN_Z should be the same size as P, or a scalar. \n');
assert(all(abs(sign_z(:)) == 1 | sign_z(:) == 0), ...
       [mfilename ':UnrecognisedSignInput'], ...
       'SIGN_Z should contain +1s, -1s or zeros. \n');

z = sign_z .* norminv(1.0 - p./2.0 - eps, 0, 1);
% [EOF]