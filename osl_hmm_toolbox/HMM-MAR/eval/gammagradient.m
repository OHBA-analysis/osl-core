function [ms,ss] = gammagradient(gamma)
% returns the average change, the highest the less stable the system is
%
% INPUTS 
% gamma:  state time courses, T X K 
%
% OUTPUTS   
% ms: mean change
% ss: sd of the change
%
% Author: Diego Vidaurre, OHBA, University of Oxford

[T, K] = size(gamma);
s = sum(abs(gamma(2:T,:)-gamma(1:T-1,:)),2);
ms = mean(s);
ss = std(s); 
