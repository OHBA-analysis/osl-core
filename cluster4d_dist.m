function stats = cluster4d_dist(S);

% subfunction of cluster4d_batch
% takes the output of cluster4d_batch and generates a null distribution
%   against which to compare cluster sizes from the true T-statistic, and
%   produces stats images
%
% input arguments:
% try, dirname = S.dirname; catch, error('Must specify S.dirname'); end
% try, tp = S.tp; catch, error('Must specify S.tp'); end
% try, np = S.np; catch, np = 1:5000; end
% try, distfile_save = S.distfile_save; catch, distfile_save = [];end
% try, distfile_load = S.distfile_load; catch, distfile_load = [];end
% try, save_images = S.save_images; catch, save_images =0; end
%
% Laurence Hunt and Tom Nichols and Mark Woolrich, 2010/11 (beta version!)


try, dirname = S.dirname; catch, error('Must specify S.dirname'); end
try, subdirname = S.subdirname; catch, error('Must specify S.subdirname'); end
try, tp = S.tp; catch, error('Must specify S.tp'); end
try, times = S.times; catch, error('Must specify S.times'); end
try, np = S.np; catch, np = 1:5000; end
try, distfile_save = S.distfile_save; catch, distfile_save = [];end
try, distfile_load = S.distfile_load; catch, distfile_load = [];end
try, save_images = S.save_images; catch, save_images =1; end
try, gridstep = S.gridstep; catch, gridstep=2;  end

nT =length(tp);

%% make clusters from true T-stat image
first = 1;
%load in the true tstat

tstats=[];
for i = 1:length(tp);
    fname = sprintf('%s/%s/%04.0f/stats%04.0f_clustere_tstat1.nii.gz',dirname,subdirname,tp(i),tp(i));
    f = nii.load(fname);
    fname_rawt = sprintf('%s/%s/%04.0f/stats%04.0f_tstat1.nii.gz',dirname,subdirname,tp(i),tp(i));
    tmp=nii.load(fname_rawt);
    if(isempty(tstats)),
       tstats=zeros([size(tmp) length(tp)]); 
    end;
    
    tstats(:,:,:,i) = tmp;

    if first %first image - read in dimensions
        Vdim = size(f);
        nV = prod(Vdim);
        X = sparse(nV,nT);
        first = 0;
    end
    X(:,i) = f(:);
end
X = clusterX(X);
Creal = unique(X);
Creal(Creal==0) = [];
if ~isempty(Creal)
    for i = 1:length(Creal)
        nVreal(i) = sum(X(find(X))==Creal(i)); %number of voxels in each indexed cluster
    end
else
    nVreal = 0;
end
clustimg = (reshape(full(X),[Vdim nT])); %cluster image for true T-stat

%% create null distribution
if isempty(distfile_load) %no distribution to load, so we create it from scratch
    dist = zeros(length(np),1);
    for p = 1:length(np) %loop over permutations
        
        first = 1;
        %% first, load in files
        for i = 1:length(tp);
            fname = sprintf('%s/%s/%04.0f/cindex%05.0f',dirname,subdirname,tp(i),np(p));
            f = nii.load(fname);
            if first %first image - read in dimensions
                nV = prod(size(f));
                X = sparse(nV,nT);
                first = 0;
            end
            X(:,i) = f(:);
        end
        
        X = clusterX(X);
        %now get the maximum number of voxels in a single cluster, and add to dist
        tmp = unique(X);
        tmp(tmp==0) = [];
        if ~isempty(tmp)
            for i = 1:length(tmp)
                nVtmp(i) = sum(X(find(X))==tmp(i)); %number of voxels in each indexed cluster
            end
        else
            nVtmp = 0;
        end
        
        dist(p) = max(nVtmp);
        clear tmp;
        fprintf('Permutation %0.0f complete...\n',p);
    end
else
    if ~iscell(distfile_load)
        error('S.distfile_load must be a cell array');
        %error('Not yet implemented loading of distributions');
    else
        for i = 1:length(distfile_load)
            dist = [];
            tmp = load(distfile_load{i});
            dist = [dist tmp.dist];
        end
    end
end

%% see where clusters lie in distribution, assign p-values and write images

if save_images

    nC = length(nVreal);
    pVreal = zeros(nC,1);
    pVimg = clustimg;
    for i = 1:nC %loop over real clusters
        pVreal(i) = mean(nVreal(i)>=dist);
        pVimg(clustimg==full(Creal(i))) = full(pVreal(i));
    end

    clustimg_fname = [dirname '/clust4d.nii.gz'];
    pVimg_fname = [dirname '/clust4d_corrp.nii.gz'];
    
    tres=1;
    nii.save(clustimg,[gridstep gridstep gridstep tres],[],clustimg_fname);
    nii.save(pVimg,[gridstep gridstep gridstep tres],[],pVimg_fname);
    
end

if ~isempty(distfile_save)
    save(distfile_save,'nVreal','clustimg','Creal','dist','times','tstats');
end

function X=clusterX(X);
nV = size(X,1);
nT = size(X,2);
% X = nV by nT matrix - input matrix & directly modified to create output matrix
% On input, X is a cluster index image *per time*; i.e. each X(:,t) gives cluster labels,
% like as output from cluster --oindex option; in particular, we do not require that cluster
% indicies are unique over time on input.
% On output, X is a cluster index image in 4D.
nci = max(X(:,1))+1; % Next cluster index  (*not* number of clusters)
for t = 2:nT
    I=X(:,t)>0; % mask of clusters for time t
    for c=unique(X(I,t))'; % cluster indicies for time t
        idx = find(X(:,t)==c);   % Voxels in cluster c at time t
        Xb  = X(idx,t-1);        % Same voxels, but cluster indicies at time t-1
        uXb = unique(Xb(Xb>0)); % Unique labels in time t-1; assume that uXb is sorted
        if length(uXb)==0
            % No relabling to do, as there is no connected cluster at time t-1
            
            % Label new cluster with unique index
            X(idx,t) = nci;
            nci = nci+1;
        else
            if length(uXb)==1
                % Exactly one cluster corresponds: the time t cluster takes on t-1's index.
                % Tag it negative, to indicate this cluster in time t is newly labeled (to avoid
                % conflict between cluster indcies in times 1:t-1 and time t).
                X(idx,t) = -uXb;
            else
                % Multiple clusters corresponds, so we need to merge
                ci = uXb(end); % use largest index for whole merged cluster
                
                % Relabel present cluster at time t, tagged negative
                X(idx,t) = -ci;
                
                % Find other voxels that need to be relabeled
                for bc = uXb(1:end-1)'
                    idx0=find(X(I,t)==-bc);      % look for previously labeled clusters in time t
                    X(I(idx0),t) = -ci;         % Relabel, but also tag it negative
                    
                    for tn = 1:t-1
                        idx0=find(X(:,tn)==bc);  % Find that cluster in time 1:t-1
                        X(idx0,tn) = ci;           % Relabel (no possibility of conflict in indicies, as 1:t-1
                        % labels are coherent, so need to flag negative)
                    end
                end
            end
        end
    end
    % clear the flags; X(:,1:t) now has consistent 4D cluster indicies
    X(:,t) = abs(X(:,t));
end


