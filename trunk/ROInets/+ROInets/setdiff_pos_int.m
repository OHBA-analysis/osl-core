function C = setdiff_pos_int(A, B)
%SETDIFF_POS_INT Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = SETDIFF_POS_INT(A, B)
%   C = A \ B   { things in A that are not in B }

%   Original by Kevin Murphy, modified by Leon Peshkin

%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 214 $
%	$LastChangedDate: 2014-07-24 12:40:42 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 14-Apr-2014 11:52:47

if isempty(A)
    C = [];
    return;
    
elseif isempty(B)
    C = A;
    return; 
    
else % both non-empty
    bits    = zeros(1, max(max(A), max(B)));
    bits(A) = 1;
    bits(B) = 0;
    
    C = A(logical(bits(A)));
end%if
end%setdiff_pos_int
% [EOF]