function dat = AB_get_source_timecourses_recon(recon,voxind)
% Giles Colclough amended L37: dat{vox,kk}(ori,timeinds,:) -> dat{vox,kk}(ori,timeinds,tri) 
if nargin == 1
  voxind = 1:length(recon.weights);
end

Ntrials = size(recon.sensor_data,3);
Ntpts   = size(recon.sensor_data,2);
Nvoxels = length(voxind);
Nori    = size(recon.weights{1},1);
NK      = numel(recon.class_samples_inds);

dat = {nan(Nori,Ntpts,Ntrials)};
dat = repmat(dat,Nvoxels,NK);

if Nvoxels > 1, ft_progress('init', 'etf'); end

for kk = 1:NK
  
  for tri = 1:Ntrials
    
    % Get subset of data for this this class and trial:
    timeinds = find(recon.class_samples_inds{kk}(1,:, tri)); % time indices for class kk
    if Ntrials == 1 && numel(recon.class_samples_inds) == 1
      sensor_data_subset = recon.sensor_data;
    else
      sensor_data_subset = recon.sensor_data(:,timeinds,tri);
    end
    
    for vox = 1:Nvoxels
      
      if Nvoxels > 1, ft_progress( ((kk-1)*(tri-1)*Nvoxels + vox) / (Nvoxels*Ntrials*NK) ); end
      
      if sum(isnan(squash(recon.weights{kk})))==0,
        if ~isempty(timeinds)
           % don't vectorise this next bit. Doesn't seem to help. GC 2014
           % (tried bsxfun, and repmat using the lightspeed mex function)
          for ori = 1:Nori
            dat{vox,kk}(ori,timeinds,tri) = recon.weights{voxind(vox),kk}(ori,:)*sensor_data_subset / sqrt(recon.wnorm_nai{voxind(vox),kk}(ori));
          end
        end
      end     
      
    end % vox
    
  end % tri
  
end % kk

if Nvoxels > 1, ft_progress('close'); end

end
% [EOF]