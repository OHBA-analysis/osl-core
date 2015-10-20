function map = glean_regress(D,regressors,mode)
% Create spatial maps via mutliple regression of HMM or ICA time courses.
%
% D = GLEAN_REGRESS(D,regressors,mode)
%
% REQUIRED INPUTS:
%   D           - Name of an SPM12 MEEG object
%   regressors  - [samples x regressors] matrix of temporal regressors
%   mode        - Type of map to create. 'pcorr' - partial correlation
%                 'tstat' - t-statistic
% 
% OUTPUTS:
%   map         - [voxels x regressors (x frequency)] spatial map
%
% Adam Baker 2015


D = spm_eeg_load(D);

F = D.nfrequencies;
if isempty(F)
    F = 1;
end

% Remove empty regressors:
regressors = double(regressors);
reg2use = ~all(regressors==0);
regressors = regressors(:,reg2use);


nRegressors = size(regressors,2);

map = nan(D.nchannels,length(reg2use),F);

if strcmp(mode,'pcorr')
    % Make std = 1
    x = bsxfun(@rdivide,regressors,std(regressors));
else
    % Remove mean
    x = bsxfun(@minus,regressors,mean(regressors));
end

if all(isinf(x(:))) % Maybe warn instead
    x(:) = 0;
end

pinvxtx = pinv(x'*x);
dof     = D.nsamples - nRegressors;

% It's slow(ish) to read data from D channel by channel so read in as much
% as possible:
blks = memblocks(size(D),1);

for iblk = 1:size(blks,1)
    
    for f = 1:F
        
        if strcmp(D.transformtype,'TF')
            Dblk = permute(D(blks(iblk,1):blks(iblk,2),f,:,1),[1 3 2]);
        else
            Dblk = D(blks(iblk,1):blks(iblk,2),:,1);
        end
        
        for v = blks(iblk,1):blks(iblk,2)
            
            y = Dblk(blks(iblk,1):blks(iblk,2) == v,:,1)';
            
            if strcmp(mode,'pcorr')
                % Normalise
                y = (y - mean(y))./std(y);
            end
            
            beta = pinvxtx * x' * y;
            
            switch mode
                case 'pcorr'
                    map(v,reg2use,f) = beta;
                case 'tstat'
                    e = y - x*beta;
                    var_e = diag(e'*e / dof);
                    t = beta ./ diag(sqrt(var_e*pinvxtx));
                    map(v,reg2use,f) = t;
            end
            
        end
        
    end
    
end
