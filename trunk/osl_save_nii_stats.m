function [statsdir,times,cluster_times]=osl_save_nii_stats( Sin )

% [statsdir,times]=osl_save_nii_stats( Sin )
%
% Outputs nii files for the specified first or group level contrasts from an OAT analysis 
% into statsdir. Can be used for results from other ( none OAT ) analyses as well.
%
% By default nii files are output in the native and 2mm resolution. Can be used for first or group level
% stats. For group level all group level contrasts are ouput for the
% specified first level contrasts.
%
% Required inputs:
% Sin.oat: oat struct resulting from an OAT analysis
%
% Sin.stats_fname: first_level or group_level results struct file name, e.g. resulting from
% an OAT analysis
% OR
% Sin.stats: first_level or group_level results struct, e.g. resulting from
% an OAT analysis
%
% Sin.stats_dir: specify ouput directory to put nii files into
% Sin.stdbrainmask_fname: std brain mask at grid resolution of the analysis
%
% Optional inputs:
% Sin.stats_dir: specify ouput directory to put nii files into
% Sin.first_level_contrasts: list of first level contrast indices
% to indicate which first level contrasts to output
% Sin.group_level_contrasts: list of group level contrast indices
% to indicate which group level contrasts to output
% should have nii files output for them
% Sin.resamp_gridstep: resolution to also resample niis to from native
% resolution. 
%
% Example first level stats call:
%
% S2=[];
% S2.oat=oat;
% S2.stats_fname=oat.first_level.results_fnames{1};
% S2.first_level_contrasts=[1,2]; % list of first level contrasts to output
% statsdir=osl_save_nii_stats(S2);
%
% Example group level stats call:
% 
% S2=[];
% S2.oat=oat;
% S2.stats_fname=oat.group_level.results_fnames;
% S2.first_level_contrasts=[1,2]; % list of first level contrasts to output
% statsdir=osl_save_nii_stats(S2);
%
% Example none OAT call:
%
% S2=[];
% S2.stdbrainmask_fname=stdbrainmask_fname;
% S2.stats_dir=stats_dir;
% stats.cope=cope;
% stats.times=times;
% S2.stats=stats; % list of first level contrasts to output
% statsdir=osl_save_nii_stats(S2);
%
% where cope is num_voxels x num_timepoints
%
% MWW 2011

global OSLDIR;

%try, tmp=Sin.oat; catch, error('Sin.oat not specfied'); end; 
%try, tmp=Sin.stats_fname; catch, error('Sin.stats_fname not specfied'); end; 

try, resamp_gridstep=Sin.resamp_gridstep; catch, resamp_gridstep=2; end; % mm

try, masksdir=[OSLDIR '/std_masks' ]; catch, error('OSLDIR not set. Run osl_startup.'); end;

try, resolve_erf_sign_ambiguity=Sin.resolve_erf_sign_ambiguity; catch, Sin.resolve_erf_sign_ambiguity=0; end;

try, tmp=Sin.freq_bin; catch, Sin.freq_bin=1; end;

have_stdcope=1;

cluster_times=[];

if(isfield(Sin,'oat')),
    if(isfield(Sin,'stats_fname')),
        Sin.stats=oat_load_results(Sin.oat,Sin.stats_fname);            
    end;

    try, statsdir=Sin.stats_dir; catch, 
        statsdir=[Sin.oat.source_recon.dirname '/' Sin.stats.fname '_dir']; 
    end; %full path of directory where stats will be stored
   
    if Sin.stats.level==1,
        stdbrainmask_fname=[Sin.oat.source_recon.dirname '/' Sin.oat.first_level.name '_mask'];
    elseif Sin.stats.level==2,
        stdbrainmask_fname=[Sin.oat.source_recon.dirname '/' Sin.oat.first_level.name '_' Sin.oat.subject_level.name '_' Sin.oat.group_level.name '_mask'];
    else
        error('Invalid stats.level');
    end;
else

    try, statsdir=Sin.stats_dir; catch, error('Sin.stats_dir not specified'); end; %full path of directory where stats will be stored
    try, stdbrainmask_fname=Sin.stdbrainmask_fname; catch, error('Sin.stdbrainmask_fname not specified'); end;
    
    if(~isfield(Sin.stats,'level')),
        Sin.stats.level=1;
    end;
    
    if(~isfield(Sin.stats,'stdcope')),
        have_stdcope=0;
    end;
end;

try, first_level_contrasts=Sin.first_level_contrasts; catch, first_level_contrasts=1:size(Sin.stats.cope,3); end; 
try, group_level_contrasts=Sin.group_level_contrasts; catch, group_level_contrasts=1:size(Sin.stats.cope,5); end; 

if(isfield(Sin.stats,'clusterstats')),
    try,corrp_thresh=Sin.cluster_stats_fpr;catch, corrp_thresh=0.05;end;
end;

tim=Sin.stats.times;

if(length(tim)>1)
    tres=tim(2)-tim(1);
else
    tres=1;
end;

if(Sin.resolve_erf_sign_ambiguity),
    Sin.stats=osl_resolve_sign_ambiguity(Sin);
    count=Sin.stats.count;
else
    count=0;
end;

if(isfield(Sin,'time_range')),
    time_indices=find(tim>=Sin.time_range(1) & tim<=Sin.time_range(2));
else
    time_indices=1:length(tim);   
end;

tim=tim(time_indices);

warning off;mkdir(statsdir);warning on;

[ mni_res ] = get_nii_spatial_res( stdbrainmask_fname );
gridstep=mni_res(1);

if(Sin.stats.level==1),

    if(size(Sin.stats.cope,4)>1)
        Sin.stats.cope=Sin.stats.cope(:,:,:,Sin.freq_bin,:);
        Sin.stats.stdcope=Sin.stats.stdcope(:,:,:,Sin.freq_bin,:);
    end;

    stdbrainmask=read_avw(stdbrainmask_fname);

    if(isfield(Sin.stats,'pseudo_zstat_var')),
        fname=[statsdir '/pseudo_zstat_var'];
        save_avw(matrix2vols(Sin.stats.pseudo_zstat_var,stdbrainmask),[fname '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, 1]);    
        %if(resamp_gridstep ~= gridstep)
        osl_resample_nii([fname '_' num2str(gridstep) 'mm'], [fname '_' num2str(resamp_gridstep) 'mm'], resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
        %end;
    end;
    
    if(isfield(Sin.stats,'sqrt_wnorm_nai')),
        fname=[statsdir '/sqrt_wnorm_nai'];
        save_avw(matrix2vols(Sin.stats.sqrt_wnorm_nai,stdbrainmask),[fname '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, 1]);    
        %if(resamp_gridstep ~= gridstep)
        osl_resample_nii([fname '_' num2str(gridstep) 'mm'], [fname '_' num2str(resamp_gridstep) 'mm'], resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
        %end;
    end;
    
    for coni=1:length(first_level_contrasts),

        con=first_level_contrasts(coni);
   
        disp(['Outputting nii files for first level contrast ' num2str(con)]);    
           
        % output nifti file for the contrast of parameter estimates
        % (COPEs) for each contrast
        fnamec=[statsdir '/cope' num2str(con)];             
        save_avw(matrix2vols(permute(Sin.stats.cope(:,time_indices,con,:),[1 2 4 3]),stdbrainmask),[fnamec '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);    
        %if(resamp_gridstep ~= gridstep)
        osl_resample_nii([fnamec '_' num2str(gridstep) 'mm'], [fnamec '_' num2str(resamp_gridstep) 'mm'], resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
        %end;

        if(have_stdcope)      
            fnamet=[statsdir '/tstat' num2str(con)];
            ts=Sin.stats.cope(:,time_indices,con,:)./Sin.stats.stdcope(:,time_indices,con,:);
            save_avw(matrix2vols(permute(ts,[1 2 4 3]),stdbrainmask),[fnamet '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);               
            %if(resamp_gridstep ~= gridstep)
            osl_resample_nii([fnamet '_' num2str(gridstep) 'mm'],[fnamet '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
            %end;
                   
            fnamev=[statsdir '/varcope' num2str(con)];                   
            save_avw(matrix2vols(permute(Sin.stats.stdcope(:,time_indices,con,:).^2,[1 2 4 3]),stdbrainmask),[fnamev '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);               
            %if(resamp_gridstep ~= gridstep)
            osl_resample_nii([fnamev '_' num2str(gridstep) 'mm'],[fnamev '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
            %end;
            
            if size(ts,2)>1,
                
                fnametsum=[statsdir '/tstat' num2str(con) '_mip'];
                time_indices_sum=find(tim>=0);
                time_indices_sum=intersect(time_indices_sum,time_indices);
                ts=abs(Sin.stats.cope(:,time_indices_sum,con,:)./Sin.stats.stdcope(:,time_indices_sum,con,:));
               
                if(0)
                    [ind]=find(ts>2.3);
                    ts2=zeros(size(ts));
                    ts2(ind)=ts(ind);                
                    ts2=sum(ts2,2)./sum(ts2>0,2);
                    ts2(isnan(ts2))=0;
                    ts=ts2;
                else
                    % output max
                    ts=max(ts,[],2);
                end;
                
                save_avw(matrix2vols(permute(ts,[1 2 4 3]),stdbrainmask),[fnametsum '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);               
                %if(resamp_gridstep ~= gridstep)
                osl_resample_nii([fnametsum '_' num2str(gridstep) 'mm'],[fnametsum '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);

            end;
        end;
        
        disp(['E.g. to view t-statistics: ']);
        disp(['!fslview ' fnamet '_' num2str(resamp_gridstep) 'mm &']);
       
    end;

elseif(Sin.stats.level==2)     
    
    if(size(Sin.stats.cope,4)>1 && ~isfield(Sin,'freq_bin'))
        error('Not implemented for multiple freq bins');
    end;
    
    if(size(Sin.stats.cope,4)>1)
        Sin.stats.cope=Sin.stats.cope(:,:,:,Sin.freq_bin,:);
        Sin.stats.stdcope=Sin.stats.stdcope(:,:,:,Sin.freq_bin,:);
    end;
    
    stdbrainmask=read_avw(stdbrainmask_fname);

    for gconi=1:length(group_level_contrasts),
        
        gcon=group_level_contrasts(gconi);
        
        for coni=1:length(first_level_contrasts),

            con=first_level_contrasts(coni);

            fnamet=[statsdir '/tstat' num2str(con) '_gc' num2str(gcon)];
            
            disp(['Outputting nii files for first level contrast ' num2str(con) ' and group contrast ' num2str(gcon)]);
            
            % output nifti files for the contrast of parameter estimates (COPEs) for each contrast
            fnamec=[statsdir '/cope' num2str(con) '_gc' num2str(gcon)]; 
            tmp=permute(Sin.stats.cope(:,time_indices,con,:,gcon),[1 2 4 3 5]);   
            save_avw(matrix2vols(tmp,stdbrainmask),[fnamec '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);    

            %if(resamp_gridstep ~= gridstep)
            osl_resample_nii([fnamec '_' num2str(gridstep) 'mm'],[fnamec '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
            %end;

            if(have_stdcope)
                tstat=permute(Sin.stats.cope(:,time_indices,con,:,gcon)./Sin.stats.stdcope(:,time_indices,con,:,gcon),[1 2 4 3 5]);           
                tstatimg=matrix2vols(tstat,stdbrainmask);
                save_avw(tstatimg,[fnamet '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);    
                if(resamp_gridstep ~= gridstep)
                    osl_resample_nii([fnamet '_' num2str(gridstep) 'mm'],[fnamet '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
                end;

                fnamev=[statsdir '/varcope' num2str(con) '_gc' num2str(gcon)];
                stat=permute(Sin.stats.stdcope(:,time_indices,con,:,gcon).^2,[1 2 4 3 5]);           
                statimg=matrix2vols(stat,stdbrainmask);
                save_avw(statimg,[fnamev '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);    
                %if(resamp_gridstep ~= gridstep)
                osl_resample_nii([fnamev '_' num2str(gridstep) 'mm'],[fnamev '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);
                %end;
                
                if size(tstat,2)>1,
                
                    fnamet_mip = [fnamet '_mip'];
                    time_indices_sum=find(tim>=0);
                    time_indices_sum=intersect(time_indices_sum,time_indices);
                    ts=abs(Sin.stats.cope(:,time_indices_sum,con,:,gcon)./Sin.stats.stdcope(:,time_indices_sum,con,:,gcon));

                    if(0)
                        [ind]=find(ts>2.3);
                        ts2=zeros(size(ts));
                        ts2(ind)=ts(ind);                
                        ts2=sum(ts2,2)./sum(ts2>0,2);
                        ts2(isnan(ts2))=0;
                        ts=ts2;
                    else
                        % output max
                        ts=max(ts,[],2);                        
                        
                        % fslmaths tstat6_full_MEGLocalSpheres_hmm0_beamform_sss0 -Tmax tstat6_full_MEGLocalSpheres_hmm0_beamform_sss0_mip
                    end;

                    stat=permute(ts,[1 2 4 3 5]);           
                    statimg=matrix2vols(stat,stdbrainmask);
                    save_avw(statimg,[fnamet_mip '_' num2str(gridstep) 'mm'],'f',[gridstep,gridstep, gridstep, tres]);    
                    %if(resamp_gridstep ~= gridstep)
                    osl_resample_nii([fnamet_mip '_' num2str(gridstep) 'mm'],[fnamet_mip '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ]);

                end;
            end;
            
            disp(['E.g. to view t-statistics: ']);
            disp(['!fslview ' fnamet '_' num2str(resamp_gridstep) 'mm &']);
       
            if(isfield(Sin.stats,'clusterstats')),
                
                disp('Outputting nii files for cluster stats');
                
                
                nC = length(Sin.stats.clusterstats{con,gcon}.Creal);
                pVreal = zeros(nC,1);
                pVimg = Sin.stats.clusterstats{con,gcon}.clustimg;
                tstat_thresh_img = zeros(size(Sin.stats.clusterstats{con,gcon}.clustimg));
                
                for i = 1:nC %loop over real clusters
                    
                    % Sin.stats.clusterstats{con,gcon}.nVreal is a list of the #voxels in each cluster 
                    % Sin.stats.clusterstats{con,gcon}.dist is the null distribution of cluster extent
                    % Sin.stats.clusterstats{con,gcon}.Creal is a list of cluster index for each cluster
                    % Sin.stats.clusterstats{con,gcon}.clustimg is a 4D image of cluster indices
                    pVreal(i) = mean(Sin.stats.clusterstats{con,gcon}.nVreal(i)>=Sin.stats.clusterstats{con,gcon}.dist); % this is 1-P-Value for each cluster
                    pVimg(Sin.stats.clusterstats{con,gcon}.clustimg==full(Sin.stats.clusterstats{con,gcon}.Creal(i))) = full(pVreal(i));
                    %if((1-pVreal(i))<corrp_thresh)
                    %    tstat_thresh_img(Sin.stats.clusterstats{con,gcon}.clustimg==full(Sin.stats.clusterstats{con,gcon}.Creal(i)))=tstatimg(Sin.stats.clusterstats{con,gcon}.clustimg==full(Sin.stats.clusterstats{con,gcon}.Creal(i)));                       
                    %end;
                    
                end;

                clustimg_fname = [fnamet '_clust4d'];
                pVimg_fname = [fnamet '_clust4d_corrp'];
                tstats_clust_fname = [fnamet '_clust4d_tstats'];
                cluster_times=Sin.stats.clusterstats{con,gcon}.times;
                
                save_avw(Sin.stats.clusterstats{con,gcon}.clustimg,[clustimg_fname '_' num2str(gridstep) 'mm'],'i',[gridstep gridstep gridstep tres]);
                save_avw(pVimg,[pVimg_fname '_' num2str(gridstep) 'mm'],'f',[gridstep gridstep gridstep tres]);
                save_avw(Sin.stats.clusterstats{con,gcon}.tstats,[tstats_clust_fname '_' num2str(gridstep) 'mm'],'f',[gridstep gridstep gridstep tres]);

                %if(resamp_gridstep ~= gridstep)
                    osl_resample_nii([clustimg_fname '_' num2str(gridstep) 'mm'],[clustimg_fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ],0);                    
                    osl_resample_nii([pVimg_fname '_' num2str(gridstep) 'mm'],[pVimg_fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'nearestneighbour',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ],0);
                    osl_resample_nii([tstats_clust_fname '_' num2str(gridstep) 'mm'],[tstats_clust_fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'nearestneighbour',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ],0);
                    runcmd(['fslmaths ' clustimg_fname '_' num2str(resamp_gridstep) 'mm -thr ' num2str(1-corrp_thresh) ' ' clustimg_fname '_thresh_' num2str(resamp_gridstep) 'mm']);
%                    runcmd(['fslmaths ' fnamet '_' num2str(resamp_gridstep) 'mm -mas ' clustimg_fname '_thresh_' num2str(resamp_gridstep) 'mm -thr ' num2str(Sin.stats.S.cluster_stats_thresh) ' ' tstat_thresh_img_fname '_' num2str(resamp_gridstep) 'mm ']);
                %end;


            end;
        end;
    end;  
    
else,
    
    error('level unsupported');
    
end

times=tim;
save([statsdir '/times'],'times','-ascii');
