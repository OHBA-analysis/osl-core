function res = bf_write_spmeeg_osl(BF, S)
% Writes out beamformer results as M/EEG dataset using different montages
% for weights normalised and unnormalised data
%
% Adam Baker
%--------------------------------------------------------------------------
if nargin == 0        
    
    none = cfg_const;
    none.tag = 'none';
    none.name = 'None';
    none.val  = {0};
    
    prefix         = cfg_entry;
    prefix.tag     = 'prefix';
    prefix.name    = 'Filename Prefix';
    prefix.help    = {'Specify the string to be prepended to the output (if relevant).'};
    prefix.strtype = 's';
    prefix.val     = {''};
    
    spmeeg      = cfg_branch;
    spmeeg.tag  = 'spmeeg_osl';
    spmeeg.name = 'SPM M/EEG dataset';
    spmeeg.val  = {prefix};
    
    res = spmeeg;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR'});

for mm  = 1:numel(modalities)

    
    % Get original channel types and units from first montage    
    montage = BF.output.montage.(modalities{mm})(1);
    chantypeorg = chantype(D, D.indchannel(montage.labelorg))';
    chanunitorg = units(D, D.indchannel(montage.labelorg))';
   
    % Apply multiple virtual montages
    for m = 1:numel(BF.output.montage.(modalities{mm}))

        montage = BF.output.montage.(modalities{mm})(m);
        montage.chantypeorg = chantypeorg;
        montage.chanunitorg = chanunitorg;

        % Online montage needs additional channel information
        for ch = 1:length(montage.labelnew)
            montage.channels(ch).label = montage.labelnew{ch};
            montage.channels(ch).type  = montage.chantypenew{ch};
            montage.channels(ch).units = montage.chanunitnew{ch};
            montage.channels(ch).bad   = 0;
        end 

        S1 = [];
        S1.montage      = montage;
        S1.keepsensors  = false;
        S1.keepothers   = false;
        S1.mode         = 'switch';


        if D.montage('getindex') ~= 0
            % If a montage is already active, then we need to chain its tra matrix
            % This functionality should not be required if using OAT because OAT only supports
            % beamforming with no online montages active. It will be encountered if directly
            % calling osl_inverse_model() after running osl_africa(). This block has been
            % tested with CTF data and could possibly throw an error with Elekta data

            fprintf('An online montage is already active, applying beamformer tra to it\n')
            old_montage = D.montage('getmontage');
            chans_to_keep = zeros(size(S1.montage.labelorg));

            for j = 1:length(chans_to_keep)
                chans_to_keep(j) = find(strcmp(old_montage.labelnew,S1.montage.labelorg{j}));
            end

            S1.montage.tra = S1.montage.tra*old_montage.tra(chans_to_keep,:);
            S1.montage.labelorg = old_montage.labelorg;
            S1.montage = rmfield(S1.montage,{'chantypeorg','chanunitorg'});
            S1.montage.name = sprintf('%s - %s',S1.montage.name,old_montage.name);
        end

        if m == 1 && mm == 1
            % Write a new MEEG object if a prefix is given
            if ~isempty(S.prefix)             
                Dnew = copy(D.montage('switch',0), fullfile(D.path,[S.prefix D.fname]));            
            else
                Dnew = D;
            end
        end;
        
        S1.D = Dnew;
        Dnew = spm_eeg_montage(S1);
        Dnew.save;

    end
end;

res.files = {fullfile(Dnew)};
