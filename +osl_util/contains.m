function ans = contains(h,q)
%
% Equivalent to: ~isempty(strfind( h, q ))
%
% Useful for older versions of Matlab.
%
% JH

    ans = ~isempty(strfind(h,q)); %#ok

end