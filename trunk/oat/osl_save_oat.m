function [oat] = osl_save_oat(oat)

% [oat] = osl_save_oat(oat)
%
% save OAT structure to the file oat.fname = [oat.source_recon.dirname '/oat_' oat.first_level.name '_' oat.subject_level.name '_' oat.group_level.name ]

if(isempty(findstr(oat.source_recon.dirname, '.oat')))
    oat.source_recon.dirname=[oat.source_recon.dirname, '.oat']; 
end;

oat.fname=[oat.source_recon.dirname '/oat_' oat.first_level.name '_' oat.subject_level.name '_' oat.group_level.name ];

save(oat.fname, 'oat');
