function [oat] = osl_load_oat(oatdir, first_level_name, subject_level_name, group_level_name)

% [oat] = osl_load_oat(oatdir, first_level_name, subject_level_name, group_level_name)
% load OAT structure with the oat.fname of [oatdir '/oat_' first_level_name '_' subject_level_name '_' group_level_name '.mat']
% OR:
%
% [oat] = osl_load_oat(oat)
%
% optional input first_level_name, subject_level_name, group_level_name
%
% MWW 2012

if nargin==0
    [oatname, oatdir]=uigetfile({'*.mat', 'MAT-files (*.mat)';'*.*', 'All Files'},'Select the OAT file');
    oat_fname = fullfile(oatdir,oatname);
else
    
    if(~isstr(oatdir)),
        
        oat=oatdir;
        oatdir=oat.source_recon.dirname;
        
        try, first_level_name=oat.first_level.name; catch first_level_name='first_level'; end;
        try, subject_level_name=oat.subject_level.name; catch subject_level_name='sub_level'; end;
        try, group_level_name=oat.group_level.name; catch group_level_name='group_level'; end;
        
    else
        
        if nargin<2,
            first_level_name='first_level';
        end;
        
        if nargin<3,
            subject_level_name='sub_level';
        end;
        
        if nargin<4,
            group_level_name='group_level';
        end;
    end;
    
    if(isempty(findstr(oatdir, '.oat')))
        oatdir=[oatdir, '.oat'];
    end;
    
    oat_fname=[oatdir '/oat_' first_level_name '_' subject_level_name '_' group_level_name '.mat'];
end
tmp=load([oat_fname]);

disp(['Loaded OAT: ' oat_fname]);

oat=tmp.oat;
oat.fname=oat_fname;

oat.source_recon.dirname=oatdir;

