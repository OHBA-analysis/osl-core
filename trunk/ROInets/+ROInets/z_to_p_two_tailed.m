function p = z_to_p_two_tailed(z)
%Z_TO_P_TWO_TAILED convert standard z-value to p-value in two-tailed test

p = 2 .* normcdf(-abs(z));
% [EOF]