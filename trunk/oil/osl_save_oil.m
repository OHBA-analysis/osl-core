function oil = osl_save_oil(oil)

% [oil] = osl_save_oil(oil)
%
% save OIL structure to the file oil.fname=[oil.source_recon.dirname '/oil_' oil.enveloping.name '_' oil.concat_subs.name '_' oil.ica.name ]

% construct name if doesn't exist
if ~isfield(oil, 'fname') || isempty(oil.fname),
    oil.fname = fullfile(oil.source_recon.dirname,        ...
                         sprintf('oil-settings-%s-%s-%s', ...
                                 oil.enveloping.name,     ...
                                 oil.concat_subs.name,    ...
                                 oil.ica.name),           ...
                          'oil.mat');
end%if

% save
ROInets.make_directory(fileparts(oil.fname));
save(oil.fname, 'oil', '-v7.3');