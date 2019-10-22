function y = osl_isactive()
%
% Tell if OSL is currently active.
% 
% JH

    y = ~isempty(getenv('OSLDIR'));

end