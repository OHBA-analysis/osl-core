function bf = tbx_cfg_bf
% Configuration file for toolbox 'Beamforming'
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tbx_cfg_bf.m 110 2014-05-15 15:06:02Z litvak.vladimir@gmail.com $

tbxdir = fileparts(mfilename('fullpath'));

if ~isdeployed, addpath(tbxdir); end

components = {
    'bf_data';
    'bf_copy';
    'bf_sources'
    'bf_features'
    'bf_inverse'
    'bf_output'
    'bf_write'
    
    };

bf = cfg_choice;
bf.tag = 'beamforming';
bf.name = 'DAiSS (beamforming)';
bf.help = {'Data analysis in source space toolbox'};

for i = 1:numel(components)
  bf.values{i} = feval(components{i});
end

% Generate the menu function automatically in case of a different directory
% name (might fail if there is no write permission)
[tbx_path, tbx_name] = fileparts(tbxdir);
if ~isequal(tbx_name, 'beamforming')
    if ~osl_util.isfile(fullfile(tbxdir, ['spm_' tbx_name '.m']))
        try
            fid = fopen(fullfile(tbxdir, ['spm_' tbx_name '.m']), 'w');
            fprintf(fid, 'function spm_%s\n\nspm_beamforming', tbx_name);
            fclose(fid);
        end
    end
end
