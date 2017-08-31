function mont = bf_output_montage(BF, S)
% Generates a montage for source extraction
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_montage.m 124 2014-11-14 13:24:42Z litvak.vladimir@gmail.com $

%--------------------------------------------------------------------------
if nargin == 0      
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for the VOI'};
    
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'MNI coordinates'; 
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the VOI in MNI coordinates'};
    pos.val = {};         
        
    radius = cfg_entry;
    radius.tag = 'radius';
    radius.name = 'Radius';
    radius.strtype = 'r';
    radius.num = [1 1];
    radius.val = {0};
    radius.help = {'Radius (in mm) for the VOIs (leave 0 for single point)'};    
    
    voidef = cfg_branch;
    voidef.tag = 'voidef';
    voidef.name = 'VOI';
    voidef.val = {label, pos, radius};
    
    vois = cfg_repeat;
    vois.tag = 'vois';
    vois.name = 'Redefine VOIs';
    vois.num  = [0 Inf];
    vois.values = {voidef};
    vois.val  = {}; 
    vois.help = {'This makes it possible to define new VOIs when the original source space was mesh or grid.',...
        'Only the sources present in the original source space can be used at this stage'}; 
   
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Summary method';
    method.labels = {'max', 'svd', 'keep'};
    method.val = {'max'};
    method.values = {'max', 'svd', 'keep'};
    method.help = {'How to summarise sources in the ROI'};
    
    mont = cfg_branch;
    mont.tag = 'montage';
    mont.name = 'Source montage';
    mont.val  = {method, vois};    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGPLANAR'});

for m  = 1:numel(modalities)    
    U        = BF.features.(modalities{m}).U; 
        
    montage          = [];
    montage.labelorg = BF.inverse.(modalities{m}).channels;
    montage.labelorg = montage.labelorg(:);
    montage.tra      = [];
    if isfield(BF.inverse.(modalities{m}), 'label')
         montage.labelnew = BF.inverse.(modalities{m}).label(:);
         montage.tra = cat(1, BF.inverse.(modalities{m}).W{:})*U';
    elseif isfield(BF.sources, 'voi') || numel(S.voidef)>0
        if isfield(BF.sources, 'voi')
            montage.labelnew = BF.sources.voi.label;
        else
            montage.labelnew = {S.voidef.label}; 
            mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI,  BF.sources.pos);
        end
        lbl = {};
        for v = 1:numel(montage.labelnew)
            if isfield(BF.sources, 'voi')
                ind = find(BF.sources.voi.pos2voi == v);
            else
                dist = sqrt(sum((mnipos-repmat(S.voidef(v).pos, size(mnipos, 1), 1)).^2, 2));
                if S.voidef(v).radius>0
                    ind = find(dist<S.voidef(v).radius);
                else
                    [minval ind] = min(dist);
                    if minval>20 % if there is nothing within 2cm something must be wrong
                        ind = [];
                    end
                end
                
                if isempty(ind)
                    error(['No sources were found close enough for VOI ' S.voidef(v).label]);
                end
            end
            
            W   = cat(1, BF.inverse.(modalities{m}).W{ind});
            
            switch S.method
                case 'max'
                    
                    Wc          = W* BF.features.(modalities{m}).C*W';  % bf estimated source covariance matrix
                    [dum, mi]   = max(diag(Wc));
                    montage.tra = [montage.tra; W(mi, :)*U'];
                    
                case 'svd'
                    %% just take top pca component for now
                    Wc          = W* BF.features.(modalities{m}).C*W'; % bf estimated source covariance matrix                                       
                       
                    [V,dum,dum]=svd(Wc);
                    montage.tra=[montage.tra;(V(:,1)'/sqrt(size(Wc, 1)))*W*U'];
                    
                case 'keep'
                    montage.tra = [montage.tra; W*U'];
                    for i = 1:size(W, 1)
                        lbl{end+1, 1} = [montage.labelnew{v} '_' num2str(i)];
                    end
            end
        end
        
        if ~isempty(lbl)
            montage.labelnew = lbl;
        end         
    else
        mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI, BF.sources.pos);
        for i = 1:size(mnipos, 1)
            w = BF.inverse.(modalities{m}).W{i};
            if ~isnan(w)
                montage.labelnew{i} = sprintf('%d_%d_%d', round(mnipos(i, :)));
                montage.tra = [montage.tra; w*U'];
            end
        end
    end
    
    montage.chantypenew = repmat({'LFP'}, length(montage.labelnew), 1);
    montage.chanunitnew = repmat({'nA*m'}, length(montage.labelnew), 1);
    
    mont.(modalities{m}) = montage;
end
