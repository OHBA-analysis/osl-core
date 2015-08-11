function out = apply_function_to_structure(f, in)
%APPLY_FUNCTION_TO_STRUCTURE applies function to every field of a structure


%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 214 $
%	$LastChangedDate: 2014-07-24 12:40:42 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 23:24:35

ff = fieldnames(in);
for iff = 1:length(ff),
    out.(ff{iff}) = f(in.(ff{iff}));
end%for

end%apply_function_to_structure