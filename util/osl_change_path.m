function [ new_full_path_fname ] = osl_change_path( old_full_path_fname, newpath )

% [ new_full_path_fname ] = osl_change_path( old_full_path_fname, newpath )

[pth fname ext]=fileparts(old_full_path_fname);
new_full_path_fname=[newpath fname ext];

end

