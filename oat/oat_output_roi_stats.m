function [stats_out,times,mni_coords_used]=oat_output_roi_stats( Sin )

% [stats_out,times]=oat_output_roi_stats( Sin )
%
% Outputs stats averaged over the specified ROI for the contrasts from an OAT analysis 
% into statsdir. Dims: ntpts x n_1st_level_con x nfreq x n_group_level_con
%
% Required inputs:
% Sin.oat: oat struct resulting from an OAT analysis
% Sin.stats_fname: first_level or group_level results struct file name, e.g. resulting from
% an OAT analysis
% Sin.mask_fname: mask to use, or set to 'all' to use all voxels
%
% Optional inputs:
% Sin.temporal_ds: temporal down sampling factor
% 
% S2=[];
% S2.oat=oat;
% S2.stats_fname=oat.first_level.results_fnames{1};
% S2.mask_fname='mask'; 
% [stats_out,times]=oat_output_roi_stats(S2);
%
% Example group level stats call:
% 
% S2=[];
% S2.oat=oat;
% S2.stats_fname=oat.group_level.results_fnames;
% S2.mask_fname='mask'; 
% [stats_out,times]=oat_output_roi_stats(S2);
%
% MWW 2011

OSLDIR = getenv('OSLDIR');

try, tmp=Sin.oat; catch, error('Sin.oat not specfied'); end; 
try, tmp=Sin.stats_fname; catch, error('Sin.stats_fname not specfied'); end; 
try, tmp=Sin.temporal_ds; catch, Sin.temporal_ds=1; end; % factor to downsample in time

Sin.stats=oat_load_results(Sin.oat,Sin.stats_fname);

try, statsdir=Sin.stats_dir; catch, statsdir=[Sin.oat.source_recon.dirname '/' Sin.stats.fname '_dir']; end; %full path of directory where stats will be stored

try, resolve_erf_sign_ambiguity=Sin.resolve_erf_sign_ambiguity; catch, Sin.resolve_erf_sign_ambiguity=0; end;

tim=Sin.stats.times;

time_indices=Sin.temporal_ds:Sin.temporal_ds:length(tim);
tim=tim(time_indices);

if(length(tim)>1)
    tres=tim(2)-tim(1);
else
    tres=1;
end;
    
stats_out=[];

if isfield(Sin,'mask_fname'),
    
    if strcmp(Sin.stats.recon_method,'none'),
        error('Must be in source space to use mask_fname');
    end;

    try,   
        mask_fname=Sin.mask_fname;
        [ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname);
    catch, 
        try,
            mni_coords=Sin.mni_coords;
        catch,
            error('A valid Sin.mask_fname or Sin.mni_coords needs to be specfied'); 
        end;
    end; 

    % find the nearest index of the beamformed voxels to the specified
    % mni_coord
    clear vox_coord
    for i=1:size(mni_coords,1),
        dists=sqrt(sum((Sin.stats.mni_coords-repmat(mni_coords(i,:)',1,size(Sin.stats.mni_coords,1))').^2,2));
        [dist,vox_coord(i,:)]=min(dists); % vox_coord is the voxel index of the beamformed voxels
    end;
    vox_coord=unique(vox_coord);
    
    mni_coords_used=Sin.stats.mni_coords(vox_coord,:);

end;

%stats_out.cope=permute(mean(Sin.stats.cope(vox_coord,:,:,:,:,:),1),[1 2 3 4 5 6]);
%stats_out.stdcope=permute(sqrt(mean(Sin.stats.stdcope(vox_coord,:,:,:,:,:).^2,1)),[1 2 3 4 5 6]);

if(Sin.resolve_erf_sign_ambiguity)
    if(size(Sin.stats.cope,4)>1),
        error('resolve_erf_sign_ambiguity not compatible with TF data');
    end;

    for c=1:size(Sin.stats.stdcope,3),
        for gc=1:size(Sin.stats.stdcope,5),
            tcs=Sin.stats.cope(vox_coord,:,c,1,gc);
            [a tc_pca]= pca(tcs,1);

            if mean(tc_pca)<0,
                tc_pca=-tc_pca;
            end;

            for ii=1:length(vox_coord),
                ccs=corrcoef(tcs(ii,:),tc_pca);
                if(ccs(1,2)<0),
                    tcs(ii,:)=-tcs(ii,:);
                end;
            end;    
            Sin.stats.cope(vox_coord,:,c,1,gc)=tcs;
            
        end;
        if isfield(Sin.stats,'lower_level_copes'),
            Sin.stats.stats.lower_level_copes{c}(vox_coord,:,:)=Sin.stats.lower_level_copes{c}(vox_coord,:,:);
        end;
    end;
end;

do_fixed_effects=0;
do_space_average=1;
dims=size(Sin.stats.cope);
dims(1)=1;
statsnew=Sin.stats;
statsnew.cope=zeros(dims);
statsnew.stdcope=zeros(dims);

for c=1:size(Sin.stats.stdcope,3),
for gc=1:size(Sin.stats.stdcope,5),

    for t=1:size(Sin.stats.stdcope,2),
    for f=1:size(Sin.stats.stdcope,4),
        S=Sin.stats.stdcope(vox_coord,t,c,f,gc).^2;
        if(do_fixed_effects)
            z=ones(length(S),1);
            y=Sin.stats.cope(vox_coord,t,c,f,gc);
            [gam, ~, covgam] = flame1(y,z,(S),do_fixed_effects);
            %[gam, covgam] = fixed_effects(y,z,diag(S));
            statsnew.cope(1,t,c,f,gc)=gam;
            statsnew.stdcope(1,t,c,f,gc)=sqrt(covgam);
        else
            if ~do_space_average,
                statsnew.cope(1,t,c,f,gc)=(Sin.stats.cope(vox_coord,t,c,f,gc));
                statsnew.stdcope(1,t,c,f,gc)=sqrt(Sin.stats.stdcope(vox_coord,t,c,f,gc).^2);
            else
                statsnew.cope(1,t,c,f,gc)=mean(Sin.stats.cope(vox_coord,t,c,f,gc),1);
                statsnew.stdcope(1,t,c,f,gc)=sqrt(mean(Sin.stats.stdcope(vox_coord,t,c,f,gc).^2,1));
            end;
        end;
        
    end;
    end;

end;
end;

statsnew.mni_coord=round(mean(mni_coords));

Sin.stats=statsnew;

Sin.stats.tstat=Sin.stats.cope./Sin.stats.stdcope;
Sin.stats.times=tim;

times=tim;

stats_out=Sin.stats;
