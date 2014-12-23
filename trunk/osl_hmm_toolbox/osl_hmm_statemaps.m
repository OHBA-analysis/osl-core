function statemaps = osl_hmm_statemaps(hmm,voxeldata,use_abs,mode)
% Computes spatial maps of state specific activity by reprojecting
% observation model variance or by fitting the HMM statepath as a regressor
% on the voxelwise data.
%
% statemaps = osl_hmm_statemaps(hmm,voxeldata,use_abs,mode)
%
% INPUT
% hmm       - the inferred HMM model structure
% voxeldata - the original full rank data or PCA components
% use_abs   - compute abs(voxeldata) before computing spatial maps
% mode      - type of spatial map to compute from the following options:
%             'var'   - outputs the variance in each state
%             'cov'   - outputs the full covariance in each state
%             'cope'  - contrast of within-state vs outside-of-state
%             'tstat' - t-statistic of within-state vs outside-of-state
%             'corr'  - correlation of data with state time course
%             'pcorr' - partial correlation of data with state time course
%
% INPUT
% statemaps - [Nvoxels x Nstates] spatial maps 
%                or 
%             [Nvoxels x Nvoxels x Nstates] matrices (mode = "cov")
%
% AB 2013



if ~exist('use_abs','var') || isempty(use_abs)
  use_abs = 0;
end

if ~exist('mode','var')
  mode = 'pcorr';
end

if exist('voxeldata','var') && ~any(strcmp(mode,{'var','cov'}))
    
    if size(voxeldata,2) ~= length(hmm.statepath)
        voxeldata = voxeldata';
    end
    
    if size(voxeldata,1) == size(hmm.MixingMatrix,1)
        Nvoxels = size(hmm.MixingMatrix,2);
    else
        Nvoxels = size(voxeldata,1);
    end
    
    statemaps = zeros(Nvoxels,hmm.K);
    
elseif strcmp(mode,'var')
    
    if isfield(hmm,'MixingMatrix')
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
    
elseif strcmp(mode,'cov')
    
    if isfield(hmm,'MixingMatrix')
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
    
    return
    
    
else
    
    error('Must specify voxeldata if not using "var" mode')
    
end

    

% Regress Viterbi path onto wholebrain results
con = cell(1,hmm.K);
for k=1:hmm.K
    con{k}=(1:hmm.K==k)';
end

cope    = zeros(Nvoxels,hmm.K);
varcope = zeros(Nvoxels,hmm.K);
c       = zeros(Nvoxels,hmm.K); 

x = zeros(length(hmm.statepath),hmm.K);

for k = 1:hmm.K
  x(:,k) = double(hmm.statepath == k);
end

if strcmp(mode,'pcorr')
  x = devar(x,1);
else
  x = demean(x,1);
end


pinvxtx = pinv(x'*x);
pinvx = pinv(x);
  


for v = 1:Nvoxels
  
    if size(voxeldata,1) == size(hmm.MixingMatrix,1)
        vdata = hmm.MixingMatrix(:,v)'*voxeldata;
    else
        vdata = voxeldata(v,:);
    end
    
  if use_abs
    if strcmp(mode,'pcorr')
      y = normalise(abs(hilbert(vdata))');
    else
      y = demean(abs(hilbert(vdata))');
    end
  else
    if strcmp(mode,'pcorr')
      y = normalise(vdata)';
    else
      y = demean(vdata)';
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
