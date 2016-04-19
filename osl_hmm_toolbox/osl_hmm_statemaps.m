function statemaps = osl_hmm_statemaps(hmm,voxeldata,use_abs,mode,assignment,de_mean,diff_contrast)
% Computes spatial maps of state specific activity by reprojecting
% observation model variance or by fitting the HMM statepath as a regressor
% on the voxelwise data.
%
% statemaps = osl_hmm_statemaps(hmm,voxeldata,use_abs,mode)
%
% INPUT
% hmm        - the inferred HMM model structure
% voxeldata  - the original full rank data or PCA components (in which case hmm
% need to contain hmm.MixingMatrix, and set data_type to 'pca')
% use_abs    - compute abs(voxeldata) before computing spatial maps
% mode       - type of spatial map to compute from the following options:
%             'var'   - outputs the variance in each state
%             'cov'   - outputs the full covariance in each state
%             'cope'  - contrast of within-state vs outside-of-state
%             'tstat' - t-statistic of within-state vs outside-of-state
%             'corr'  - correlation of data with state time course
%             'pcorr' - partial correlation of data with state time course
%             'conn'  - compute the connectivity profile of each state
% data_type  - type of data in voxeldata
%            - 'voxel'
%            - 'pca'
%
% assignment - 'hard' - use hard state assignment (Viterbi path, default)
%              'soft' - use probabilistic state assignment
% INPUT
% statemaps - [Nvoxels x Nstates] spatial maps 
%                or 
%             [Nvoxels x Nvoxels x Nstates] matrices (mode = "cov")
%
% AB 2013

if ~exist('diff_contrast','var')
    diff_contrast=1;
end

if ~exist('de_mean','var')
    de_mean = 1; 
end;

if ~exist('assignment','var')
    assignment = 'hard'; 
end;

if ~exist('use_abs','var') || isempty(use_abs)
  use_abs = 0;
end

if ~exist('mode','var')
  mode = 'pcorr';
end

if ~exist('data_type')
    data_type='voxel';
end

if strcmp(data_type,'pca')
    apply_mixing_matrix=1;
elseif strcmp(data_type,'voxel')
    apply_mixing_matrix=0;
else
    error('Illegal datatype');
end

if exist('voxeldata','var') && ~any(strcmp(mode,{'var','cov','conn'}))
    
    if size(voxeldata,2) ~= length(hmm.statepath)
        voxeldata = voxeldata';
    end
    
    if apply_mixing_matrix && size(voxeldata,1) == size(hmm.MixingMatrix,1)
        Nvoxels = size(hmm.MixingMatrix,2);
    else
        Nvoxels = size(voxeldata,1);
    end
    
    statemaps = zeros(Nvoxels,hmm.K);
    
elseif strcmp(mode,'var')
    
    if apply_mixing_matrix && isfield(hmm,'MixingMatrix')
        statemaps = zeros(size(hmm.MixingMatrix,2),hmm.K);
        for k = 1:hmm.K
            statemaps(:,k) = diag(hmm.MixingMatrix'*hmm.state(k).Cov*hmm.MixingMatrix);
        end
    else
        statemaps = zeros(size(hmm.state(1).Cov,1),hmm.K);
        for k = 1:hmm.K
            statemaps(:,k) = diag(hmm.state(k).Cov);
        end
    end
    
    return
    
elseif any(strcmp(mode,{'cov','conn'}))
    
    if exist('voxeldata','var') && ~isempty(voxeldata)
        statemaps = zeros(size(voxeldata,1),size(voxeldata,1),hmm.K);
        for k = 1:hmm.K
            statemaps(:,:,k) = cov(voxeldata(:,hmm.statepath==k)');
        end
    else
        if apply_mixing_matrix && isfield(hmm,'MixingMatrix')
            statemaps = zeros(size(hmm.MixingMatrix,2),size(hmm.MixingMatrix,2),hmm.K);
            for k = 1:hmm.K
                statemaps(:,:,k) = hmm.MixingMatrix'*hmm.state(k).Cov*hmm.MixingMatrix;
            end
        else
            statemaps = zeros(size(hmm.state(1).Cov,1),size(hmm.state(1).Cov,1),hmm.K);
            for k = 1:hmm.K
                statemaps(:,:,k) = hmm.state(k).Cov;
            end
        end
        
    end
    
    if strcmp(mode,'conn')
        
        Nvoxels = size(statemaps,1);
        Conn = zeros(Nvoxels,hmm.K);
        
        for k = 1:hmm.K
            statemaps(:,:,k) = corrcov(statemaps(:,:,k));
        end
    
        statemaps(repmat(logical(eye(Nvoxels)),[1,1,hmm.K])) = nan;
        
        for vox = 1:size(statemaps,1)
            for k = 1:hmm.K
                centroid = mean(statemaps(vox,1:Nvoxels~=vox,1:hmm.K~=k),3);
                cp = statemaps(vox,1:Nvoxels~=vox,k);
                Conn(vox,k) = norm((centroid - cp),2); 
            end
        end
        
        statemaps = Conn;

    end
    

    return
    
else
    
    error('Must specify voxeldata if not using "var" mode')
    
end


% Regress Viterbi path onto wholebrain results
con = cell(1,hmm.K);

if strcmp(mode,'cope') | strcmp(mode,'tstat')
    for k=1:hmm.K
        if diff_contrast
            con{k}=-1*ones(hmm.K,1)*(1/(hmm.K-1));
            con{k}(k)=1;
        else
            con{k}=(1:hmm.K==k)';
        end
    end
else
    for k=1:hmm.K
        con{k}=(1:hmm.K==k)';
    end
end

cope    = zeros(Nvoxels,hmm.K);
varcope = zeros(Nvoxels,hmm.K);
c       = zeros(Nvoxels,hmm.K); 

x = zeros(length(hmm.statepath),hmm.K);

switch assignment
    case 'hard'
        for k = 1:hmm.K
            x(:,k) = double(hmm.statepath == k);
        end
    case 'soft' 
        x = hmm_sub.statepath_soft;
end

if strcmp(mode,'pcorr')
  x = devar(x,1);
  x(isnan(x))=0;
else
  if de_mean
    x = demean(x,1);
  end
end


pinvxtx = pinv(x'*x);
pinvx = pinv(x);
  


for v = 1:Nvoxels
  
    if apply_mixing_matrix && size(voxeldata,1) == size(hmm.MixingMatrix,1)
        vdata = hmm.MixingMatrix(:,v)'*voxeldata;
    else
        vdata = voxeldata(v,:);
    end
    
  if use_abs
    if strcmp(mode,'pcorr')
      y = normalise(abs(hilbert(vdata))');
    else
      if de_mean
        y = demean(abs(hilbert(vdata))');
      else
        y= abs(hilbert(vdata))';
      end
    end
  else
    if strcmp(mode,'pcorr')
      y = normalise(vdata)';
    else
      if de_mean
        y = demean(vdata)';
      else
        y=vdata';
      end
    end
  end
    
  if strcmp(mode,'corr')
    c(v,:) = corr(x,repmat(y,1,hmm.K));
  else
    [cope(v,:),varcope(v,:)] = glm_fast_for_meg(y,x,pinvxtx,pinvx,con,0);
    tstat = cope ./ sqrt(varcope);
  end

end


switch mode
  case {'cope','pcorr'}
    statemaps = cope;
  case 'tstat'
    statemaps = tstat;
  case 'corr'
    statemaps = c;
end
