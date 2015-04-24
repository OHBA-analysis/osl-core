function [ D ] = oat_get_sensordata( oat_results, oat )

% [ D ] = oat_get_sensordata( oat_results )
% [ D ] = oat_get_sensordata( oat_results_fname, oat )
%
% get D and update path in case OAT dir has been moved
% 
% MWW 2014

if isstr(oat_results)
    if nargin <2
        error('Need to specify oat if passing in oat_results unloaded as a file name');
    end;
    
    oat_results = oat_load_results(oat,oat_results);
end;

if(isfield(oat_results,'BF'))
    D=fullfile(oat_results.BF.data.D);
else,
    D=fullfile(oat_results.D_sensor_data);
end;

[D_path, D_nme, D_ext]=fileparts(D);

oatdirname=oat_results.source_recon.dirname;

% if needs be adjust path to the oatdir name (e.g. will be
% necessary if the oatdir has been moved or renamed)
if(~strcmp([oatdirname],D_path)),  
    D=[oatdirname,'/',D_nme,D_ext];

    D=spm_eeg_load(D);
    
    D=path(D,oatdirname);
    
    D.save;
else
    D=spm_eeg_load(D);    
end;

% make sure we are using sensor space montage:
D=montage(D,'switch',0);

%D = osl_correct_meeg_paths(D);

end

