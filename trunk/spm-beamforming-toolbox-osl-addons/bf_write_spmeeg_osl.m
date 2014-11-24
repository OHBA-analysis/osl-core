function res = bf_write_spmeeg_osl(BF, S)
% Writes out beamformer results as M/EEG dataset using different montages
% for weights normalised and unnormalised data
%
% Adam Baker
%--------------------------------------------------------------------------
if nargin == 0
    
    modality         = cfg_menu;
    modality.tag     = 'modality';
    modality.name    = 'Modality';
    modality.help    = {'What modality to output'};
    modality.labels  = {
        'MEG'
        'MEGPLANAR'
        'EEG'
        }';
    modality.values  = {
        'MEG'
        'MEGPLANAR'
        'EEG'
        }';
    modality.val = {'MEG'};
    
    none = cfg_const;
    none.tag = 'none';
    none.name = 'None';
    none.val  = {0};
    
    prefix         = cfg_entry;
    prefix.tag     = 'prefix';
    prefix.name    = 'Filename Prefix';
    prefix.help    = {'Specify the string to be prepended to the output (if relevant).'};
    prefix.strtype = 's';
    prefix.num     = [1 Inf];
    prefix.val     = {'B'};
    
    spmeeg      = cfg_branch;
    spmeeg.tag  = 'spmeeg_osl';
    spmeeg.name = 'SPM M/EEG dataset';
    spmeeg.val  = {modality, prefix};
    
    res = spmeeg;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

% Get original channel types and units from first montage    
montage = BF.output.montage.(S.modality)(1);
chantypeorg = chantype(D, D.indchannel(montage.labelorg))';
chanunitorg = units(D, D.indchannel(montage.labelorg))';

% Apply multiple virtual montages
for m = 1:numel(BF.output.montage.(S.modality))
    
    montage = BF.output.montage.(S.modality)(m);

    montage.chantypeorg = chantypeorg;
    montage.chanunitorg = chanunitorg;
    
    S1 = [];
    S1.montage      = montage;
    S1.prefix       = S.prefix; % ignored for online
    S1.keepsensors  = false;
    S1.keepothers   = false;
    S1.mode         = 'switch';
    
    %if m == 1
    %  D = copy(D, [S.prefix D.fname]);
    %end
    
    S1.D = D;
    D = spm_eeg_montage(S1);
    D.save;
    
end

res.files = {fullfile(D)};
