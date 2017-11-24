function mont = bf_output_montage_osl(BF, S)
% Generates a montage for source extraction with weights normalisation
%
% Adam Baker & MWW
%--------------------------------------------------------------------------
if nargin == 0
        
    normalise        = cfg_menu;
    normalise.tag    = 'normalise';
    normalise.name   = 'Apply weights normalisation';
    normalise.labels = {'yes', 'no','keep both'};
    normalise.values = {'yes','no','both'};
    normalise.val    = {'both'};
    
    mont = cfg_branch;
    mont.tag = 'montage_osl';
    mont.name = 'Source montage with weights normalisation';
    mont.val  = {normalise}; 
    
    return
    
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR'});

for m  = 1:numel(modalities)
    
    montage   = [];
    montage.labelorg = BF.inverse.(modalities{m}).channels(:);    
    
    %%%%%%
    % get weights to work out sizes of containers
    if isfield(BF.inverse.(modalities{m}),'class'),
        % using multi-class festures (see bf_features)   
        W=BF.inverse.(modalities{m}).class{1}.W;
        nclasses=length(BF.inverse.(modalities{m}).class);
    else
        W=BF.inverse.(modalities{m}).W;
        nclasses=1;
    end;
    
    Ncoords = length(W);

    Nori = size(W{1},1);
    %%%%%
    
    %%%%%
    % Set up separate tra matrices for weights and normalised weights
    tra_w     = zeros(Ncoords*Nori,size(W{1}(1,:),2));
    tra_wnorm = zeros(Ncoords*Nori,size(W{1}(1,:),2));
    
    lbl = cell(Nori*Ncoords,1);
       
    montage.labelnew = cell(Ncoords,1);
    %%%%%
        
    for kk=1:nclasses,

        for ind = 1:Ncoords % for all points in space

            % Get weights
            if isfield(BF.inverse.(modalities{m}),'class'),
                % using multi-class festures (see bf_features)   
                W=cat(1,BF.inverse.(modalities{m}).class{kk}.W{ind});
            else
                W=cat(1,BF.inverse.(modalities{m}).W{ind});
            end;

            % set label for this coord
            montage.labelnew{ind} = mat2str(BF.sources.pos(ind,:),3);

            for ii = 1:Nori % for all orientations
                % Index into tra is a combination of voi & ori - e.g.
                % v1o1,v1o2,v2o1...
                tra_ind = (ind-1)*Nori + ii;
                lbl{tra_ind} = [montage.labelnew{ind} ', ' num2str(ii)];
                tra_w(tra_ind,:)     = W(ii,:);
                tra_wnorm(tra_ind,:) = W(ii,:)./(sqrt(W(ii,:)*W(ii,:)'));
            end;

        end
    
        if ~isempty(lbl)
            montage.labelnew = lbl;
        end

        montage.chantypenew = repmat({'LFP'}, length(montage.labelnew), 1);
        montage.chanunitnew = repmat({'nA*m'}, length(montage.labelnew), 1);

        switch S.normalise
            case {'no'}
                montage.tra = tra_w;
                montage.name = sprintf('Source space (%s) without weights normalisation, class %d',modalities{m},kk);
                mont.(modalities{m})(kk) = montage;
            case {'yes'}
                montage.tra = tra_wnorm;
                montage.name = sprintf('Source space (%s) with weights normalisation, class %d',modalities{m},kk);
                mont.(modalities{m})(kk) = montage;
            case {'both'} % Save two montages, weights and weights normalised
                montage.tra = tra_w;
                montage.name = sprintf('Source space (%s) without weights normalisation, class %d',modalities{m},kk);
                mont.(modalities{m})(kk) = montage;
                montage.tra = tra_wnorm;
                montage.name = sprintf('Source space (%s) with weights normalisation, class %d',modalities{m},kk);
                mont.(modalities{m})(nclasses+kk) = montage;    
        end
    end;

end