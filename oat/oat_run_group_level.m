function [ results_fnames current_level_results ] = osl_run_group_level( oat )

% [ results_fnames ] = osl_run_group_level( oat )
%
% takes in an OAT, which needs to be setup by calling oat=osl_setup_oat(S), struct
% and runs group level subjectwise GLM
%
% This function should normally be called using osl_run_oat(oat);
%
% MWW 2011

OSLDIR = getenv('OSLDIR');

lower_level=oat.subject_level;
current_level=oat.group_level;

store_lower_level_copes=oat.group_level.store_lower_level_copes;

try, glm_method=current_level.glm_method; catch, glm_method='ols'; end;

use_flame1=0;
if(~strcmp(glm_method,'ols'))
    use_flame1=1;
end;

if(current_level.use_tstat && ~strcmp(glm_method,'ols')),
    current_level.use_tstat=0;
    warning('Switching to using first level copes rather than t-stats, since a flame1 style analysis is being used');
end;

if(~strcmp(glm_method,'fixed_effects') && length(oat.first_level.results_fnames)<=size(oat.group_level.group_design_matrix,1)+2)
    error('Not enough subjects for design matrix to be fitted using mixed effects (e.g. ols).');
end;

disp(['Loading in lower level stats']);

current_level_results=[];

ft_progress('init', 'etf');

for subi=1:length(current_level.subjects_to_do),

    ft_progress(subi/length(current_level.subjects_to_do));

    sub=oat.subject_level.subjects_to_do(current_level.subjects_to_do(subi));  
    subinds(subi)=sub;
    
    % load in single subject stats
    lower_level_results=oat_load_results(oat, lower_level.results_fnames{sub});
        
    if(subi==1)
        sc=[size(lower_level_results.cope) 1 1];
        cope=zeros(sc(1),length(current_level.subjects_to_do),sc(2),sc(3),sc(4),'single');
        stdcope=inf(sc(1),length(current_level.subjects_to_do),sc(2),sc(3),sc(4),'single');       
        
        current_level_results.fname=[oat.first_level.name '_' lower_level.name '_' current_level.name];

        current_level_results.recon_method=lower_level_results.recon_method;

        if strcmp(current_level_results.recon_method,'none'),
            disp(['Using whole of sensor space']);
            current_level_results.mask_indices_in_lower_level=1:size(cope,1);
                       
            current_level_results.chanind=lower_level_results.chanind;
            current_level_results.chanlabels=lower_level_results.chanlabels;
            current_level_results.chantype=lower_level_results.chantype;
                       
        else,

            if(isfield(oat.first_level,'mni_coords'))
                
                disp('Using MNI coordinates (not mask) at group level.')
                disp('If these MNI coordinates are spatially sparse (e.g. to select ROIs) be sure to turn all spatial smoothing off at the group level');
                if isfield(oat.group_level,'mni_coords')
                    current_level_mni_coord = oat.group_level.mni_coords;
                else
                    current_level_mni_coord = oat.first_level.mni_coords;
                end
                
                lower_level_mni_coord=lower_level_results.mni_coords;
                
                [~, iA, mask_indices]=intersect(current_level_mni_coord, lower_level_mni_coord,'rows');
                
                % sort indices to be in the order of the lower level mni_coords:
                [~, gg]=sort(iA);
                indices_in_lower_level=mask_indices(gg);
                
                current_level_results.mask_indices_in_lower_level = indices_in_lower_level;
                current_level_results.mni_coords = lower_level_mni_coord(indices_in_lower_level,:);
                
                % the mni coordinates from the specified list
                str=['Desired MNI coordinates ' ];
                for vox=1:size(current_level_mni_coord,1),
                    str=[str ', [' num2str(current_level_mni_coord(vox,:)) ']'];
                end;
                disp(str);
                
                % the mni coordinates from the lower level mask
                str=['Nearest dipoles from lower level at MNI coordinates ' ];
                for vox=1:length(indices_in_lower_level),
                    str=[str ', [' num2str(lower_level_mni_coord(indices_in_lower_level(vox),:)) ']'];
                end;
                disp(str);
                               

            elseif oat.first_level.parcellation.do

                current_level_results.mask_indices_in_lower_level = 1:(size(lower_level_results.D_sensor_data,1) - 1);
                current_level_results.mni_coords = lower_level_results.D_sensor_data.parcellation.mni_coords;

                % need to also pass D through to group leve results so that
                % parcellation info is available
                current_level_results.D_sensor_data=lower_level_results.D_sensor_data;
            else
                
                current_level_results.gridstep=lower_level_results.gridstep;
                
                S=[];
                
                S.lower_level_mask_fname=[oat.source_recon.dirname filesep oat.first_level.name '_' lower_level.name '_mask.nii.gz'];
                S.current_level_mask_fname=[oat.source_recon.dirname filesep oat.first_level.name '_' lower_level.name '_' current_level.name '_mask.nii.gz'];
                
                S.current_level=current_level;
                S.lower_level_mni_coord=lower_level_results.mni_coords;
                S.lower_level_gridstep=lower_level_results.gridstep;
                
                [current_level_results.mask_indices_in_lower_level,current_level_results.mni_coords,currmask]=setup_mask_indices(S);
                
                % group level mask for doing spatial averaging
                stdbrain_fname  = S.current_level_mask_fname;
                curr_level_mask_fname = stdbrain_fname;
                lower_level_stdbrain_fname = S.lower_level_mask_fname;
                
            end
            
        end;
    end;
    
    cope(:,subi,:,:,:)=lower_level_results.cope; % nvox x nsub x ntpts x ncons x nfreqs

    stdtmp=std(squash(cope(:,subi,:,:,:)));

    stdcope(:,subi,:,:,:)=lower_level_results.stdcope; 
    
    for con=1:size(lower_level_results.cope,4),
        % check for infs or zeros in stdcope
        if(sum(squash(isinf(lower_level_results.stdcope)))>0 || any(squash(lower_level_results.stdcope)==0)),
            warning(['Subject ' num2str(sub) ' has infinite or zero stdcope for first-level contrast ' num2str(con) '. This may be due to there being zero trials for one of the first level conditions.']);
            if(strcmp(glm_method,'ols'))
                error(['Can not continue with OLS.']);
            end;
        end;
    end;
end;

ft_progress('close');


lower_level_results=rmfield(lower_level_results,'cope');
lower_level_results=rmfield(lower_level_results,'stdcope');

%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%
%% now setup to do GLM stats

num_subjects=size(current_level.group_design_matrix,2);

x=current_level.group_design_matrix';
pinvx=pinv(x);
pinvxtx=pinv(x'*x);

group_contrast=current_level.group_contrast;

if(size(cope,2)~=num_subjects)
    error('mismatched data, size(cope,2)~=num_subjects. Check your design matrix is compatible with your subject_todo list.');
end;

for gc=1:length(current_level.group_contrast),
    if(size(current_level.group_design_matrix,1)~=length(current_level.group_contrast{gc})),
        error(['Mismatched group design matrix and group contrast ' num2str(gc) '. Number of regressors are not the same.']);
    end;
end;

if(length(current_level.group_contrast_name)~=length(current_level.group_contrast)),
    error(['Mismatched number of group contrast names and group contrasts']); 
end;

ncons=size(cope,4);
nfreq=size(cope,5);

if(current_level.use_tstat)
    cope=cope./stdcope;
end;

%mincope=abs(min(squash(cope)))+1e-10;
%cope=log(cope+mincope);

%cope=log(cope);

%if(~use_flame1),
%    clear stdcope; % clear up some memory
%end;

if(size(cope,3)~=length(lower_level_results.times)),
    error('size(cope,3)~=length(lower_level_results.times)');
end;

if(size(cope,5)~=length(lower_level_results.frequencies)),
    error('size(cope,5)~=length(lower_level_results.frequencies)');
end;

%% setup time and frequency windows/indices for using later

if(isempty(current_level.time_range))
    current_level.time_range=lower_level_results.first_level.time_range;
end;

time_indices=lower_level_results.times>=current_level.time_range(1) & lower_level_results.times<=current_level.time_range(2);
current_level_results.times=lower_level_results.times(time_indices);
current_level_results.time_indices_into_lower_level_times=time_indices;
current_level_results.lower_level_times=lower_level_results.times;

if(length(lower_level_results.times)>1)
    tres=lower_level_results.times(2)-lower_level_results.times(1);
else
    tres=1;
end;

freq_indices=1;
if(lower_level_results.frequencies>0),
    if(isempty(current_level.freq_range))
        current_level.freq_range=lower_level_results.first_level.tf_freq_range;
    end;
    freq_indices=lower_level_results.frequencies>=current_level.freq_range(1) & lower_level_results.frequencies<=current_level.freq_range(2);
end;
current_level_results.frequencies=lower_level_results.frequencies(freq_indices);
current_level_results.frequency_ranges=lower_level_results.frequency_ranges(freq_indices,:);

if(current_level.time_average),
    current_level_results.times=mean(current_level_results.times);
    tres=1;
end;

if(current_level.freq_average),
    current_level_results.frequencies=mean(current_level_results.frequencies);
end;

Nvoxels_out=length(current_level_results.mask_indices_in_lower_level);
if(current_level.space_average),
    Nvoxels_out=1;
    if ~strcmp(current_level_results.recon_method,'none'),
        current_level_results.mni_coords=mean(current_level_results.mni_coords,1);
    end;
end;

ntpts=length(current_level_results.times);
nfreq=length(current_level_results.frequencies);

%% results containers
current_level_results.cope=zeros(Nvoxels_out,ntpts,ncons,nfreq,length(group_contrast),'single');
current_level_results.stdcope=inf(Nvoxels_out,ntpts,ncons,nfreq,length(group_contrast),'single'); % nvox x ntpts x ncons x nfreqs x ngcons

%% do a 1st level contrast at a time
for cc=1:length(current_level.first_level_contrasts_to_do),
    
    c=current_level.first_level_contrasts_to_do(cc);
    
    disp(['Processing Group stats for first level contrast ' num2str(c) ', ' lower_level_results.first_level_contrast_name{c}]);
    disp(['Number ' num2str(cc) ' of ' num2str(length(current_level.first_level_contrasts_to_do))]);
    
    disp(['Preparing data']);
    
    %% temporary container to be reused for each contrast
    cope_smooth_lower_level=zeros(size(cope,1),size(cope,2),length(lower_level_results.times),length(lower_level_results.frequencies),'single');% nvox x nsub x ntpts x nfreqs
    stdcope_smooth_lower_level=inf(size(cope,1),size(cope,2),length(lower_level_results.times),length(lower_level_results.frequencies),'single');% nvox x nsub x ntpts x nfreqs
    
    %% extract lower level data:
    for voxi=1:size(cope,1),
        
        for sub=1:num_subjects,
            
            cope_smooth_lower_level(voxi,sub,:,:)=permute(cope(voxi,sub,:,c,:),[1,2,3,5,4]); % nvox x nsub x ntpts x nfreqs
            stdcope_smooth_lower_level(voxi,sub,:,:)=permute(stdcope(voxi,sub,:,c,:),[1,2,3,5,4]); % nvox x nsub x ntpts x nfreqs
            
        end;
    end;
    
    %% check for subjects with extremely high lower level stdcopes            
    avg_stdcope=permute(mean(mean(mean(stdcope_smooth_lower_level,1),3),4),[2,1,3,4]);
    outlier_avg_stdcope=(((avg_stdcope-mean(avg_stdcope))/std(avg_stdcope))>3);
    if ~isempty(find(outlier_avg_stdcope)),
        for sub=1:length(outlier_avg_stdcope),
            if(outlier_avg_stdcope(sub))
                warning(['Subject ' num2str(subinds(sub)) ' has high stdcope for first-level contrast ' num2str(con) '. This may be due to there being zero trials for one of the first level conditions. Consider removing this subject from the group analysis.']);           
            end; 
        end;
    end;

    %% do spatial smoothing
    if(current_level.spatial_smooth_fwhm>0),
        
        disp(['Spatial smoothing with FWHM=' num2str(current_level.spatial_smooth_fwhm) ' mm']);
        if strcmp(current_level_results.recon_method,'none'),
            if(sub==1)
                warning('Not implemented for sensor space - no spatial smoothing applied');
            end;
        else
            
            
            for sub=1:num_subjects,
                for f=1:nfreq,
                    
                    cope_smooth_lower_level(:,sub,:,f)=smooth_vol(permute(cope_smooth_lower_level(:,sub,:,f),[1 3 2 4]), lower_level_stdbrain_fname, current_level.spatial_smooth_fwhm);
                    stdcope_smooth_lower_level(:,sub,:,f)=smooth_vol(permute(stdcope_smooth_lower_level(:,sub,:,f),[1 3 2 4]), lower_level_stdbrain_fname, current_level.spatial_smooth_fwhm);
                    
                end;
            end;
                        
        end;
    end;
    
    % vox=10; figure;plot(squeeze(cope_smooth_lower_level(vox,1,:,1)));
    
    %% do temporal smoothing
    if(current_level.time_smooth_std>0),
        disp(['Temporal smoothing with std=' num2str(current_level.time_smooth_std) ' secs']);
        
        for voxi=1:size(cope_smooth_lower_level,1),
            
            for sub=1:num_subjects,
                
                dat=permute(cope_smooth_lower_level(voxi,sub,:,:),[3,5,1,2,4]);
                datstd=permute(stdcope_smooth_lower_level(voxi,sub,:,:),[3,5,1,2,4]).^2;
                
                if(nfreq>1)
                    if(sub==1)
                        warning('Temporal smoothing (current_level.time_smooth_std>0) not implemented for time-freq analysis, so none is applied');
                    end;
                    
                    dat2=dat;
                    dat2std=datstd;
                    
                else
                    
                    f = fftshift(osl_gauss(current_level.time_smooth_std/tres,1,length(datstd))');
                    
                    dat2 = fftconv(dat,f);
                    dat2std = fftconv(datstd,f);
                    
                end;
                
                cope_smooth_lower_level(voxi,sub,:,:)=dat2;
                stdcope_smooth_lower_level(voxi,sub,:,:)=dat2std;
                
                % figure;plot(dat);ho;plot(dat2,'r');
            end;
        end;
    end;
    
    %% now do windowing for all dimensions
    cope_smooth_lower_level=cope_smooth_lower_level(current_level_results.mask_indices_in_lower_level,:,time_indices,freq_indices);
    stdcope_smooth_lower_level=stdcope_smooth_lower_level(current_level_results.mask_indices_in_lower_level,:,time_indices,freq_indices);
    
    %if(current_level.use_voxelwise_stdnorm)
    %    vars=repmat(mean(mean(stdcope_smooth_lower_level.^2,2),3),[1,size(stdcope_smooth_lower_level,2),size(stdcope_smooth_lower_level,3)]);
    %    cope_smooth_lower_level=cope_smooth_lower_level./sqrt(vars);
    %    stdcope_smooth_lower_level=stdcope_smooth_lower_level./sqrt(vars);
    %end;
    
    %% now do averaging over dimensions
    if(current_level.space_average)
        % average over space
        cope_smooth_lower_level=mean(cope_smooth_lower_level,1);
        stdcope_smooth_lower_level=sqrt(mean(stdcope_smooth_lower_level.^2,1));
    end;
    
    if(current_level.time_average)
        % average over time bins
        cope_smooth_lower_level=mean(cope_smooth_lower_level,3);
        stdcope_smooth_lower_level=sqrt(mean(stdcope_smooth_lower_level.^2,3));
        
        if(current_level.group_varcope_time_smooth_std>0),
            warning('Doing time averaging. This means that there will be no variance smoothing of the between-subject variance (as specified by group_varcope_time_smooth_std).');
        end;
        
    end;
    
    if(current_level.freq_average)
        % average over frequency bins
        cope_smooth_lower_level=mean(cope_smooth_lower_level,4);
        stdcope_smooth_lower_level=sqrt(mean(stdcope_smooth_lower_level.^2,4));
    end;
    
    if(current_level.use_tstat && ~use_flame1)
        %cope_smooth_lower_level=cope_smooth_lower_level./stdcope_smooth_lower_level;
    end;
    
    %% now do the GLM stats
    disp(['Computing group statistics']);
    ft_progress('init', 'etf');
    
    for vox=1:size(cope_smooth_lower_level,1),
        
        ft_progress(vox/size(cope_smooth_lower_level,1));
        
        for f=1:size(cope_smooth_lower_level,4),
            
            datf=permute(cope_smooth_lower_level(vox,:,:,f),[2 3 1 4 5]);
            
            %mn=mean(squash(datf(:,current_level_results.times<0)));
            mn=0;
            
            for t=1:size(cope_smooth_lower_level,3),
                
                dat=datf(:,t)-mn;
                
                % stddat=std(dat); dat=dat/stddat;
                stddat=1;
                
                if(use_flame1)
                    Sdat=permute(stdcope_smooth_lower_level(vox,:,t,f),[2 1 3 4 5]);
                    Sdat=(Sdat/stddat).^2;
                end;
                
                if(current_level.use_robust_glm)
                    if(strcmp(glm_method,'ols')),
                        [b,statsrf]=robustfit(x,squeeze(dat)','bisquare',4.685,'off'); % fits GLM robustly, Y=XB+E
                        
                        for gc=1:length(group_contrast),
                            gcon=group_contrast{gc}';
                            current_level_results.cope(vox,t,c,f,gc)=gcon*b; % num_group_contrasts x num_contrasts x num_timepoints x num_freqs x num_vox
                            current_level_results.stdcope(vox,t,c,f,gc)=sqrt(gcon*statsrf.covb*gcon');
                        end;
                    else
                        error(['glm_method ' glm_method ' is unsupported with robust GLM']);
                    end;
                else
                    if(strcmp(glm_method,'ols')),
                        [copeout, varcopeout, coapeout, dof]=glm_fast_for_meg(dat,x,pinvxtx,pinvx,group_contrast,0);
                        current_level_results.stdcope(vox,t,c,f,:)=sqrt(varcopeout);
                        current_level_results.cope(vox,t,c,f,:)=copeout;
                    elseif(strcmp(glm_method,'flame1')),
                        [b, beta, covgam]=flame1(squeeze(dat),x,squeeze(Sdat));
                        
                        for gc=1:length(group_contrast),
                            gcon=group_contrast{gc}';
                            current_level_results.cope(vox,t,c,f,gc)=gcon*b; % num_group_contrasts x num_contrasts x num_timepoints x num_freqs x num_vox
                            current_level_results.stdcope(vox,t,c,f,gc)=sqrt(gcon*covgam*gcon');
                        end;
                    elseif(strcmp(glm_method,'fixed_effects')),
                        [b, beta, covgam]=flame1(squeeze(dat),x,squeeze(Sdat),1);
                        
                        for gc=1:length(group_contrast),
                            gcon=group_contrast{gc}';
                            current_level_results.cope(vox,t,c,f,gc)=gcon*b; % num_group_contrasts x num_contrasts x num_timepoints x num_freqs x num_vox
                            current_level_results.stdcope(vox,t,c,f,gc)=sqrt(gcon*covgam*gcon');
                        end;
                    else
                        error(['glm_method ' glm_method ' unsupported']);
                    end;
                end;
                
            end;
            
            
            %% do temporal smoothing of stds
            if(current_level.group_varcope_time_smooth_std>0 && ~current_level.time_average)
                for gc=1:length(group_contrast),
                    dat=permute(current_level_results.stdcope(vox,:,c,f,gc),[2,1,3,4,5]);
                    
                    st=current_level.group_varcope_time_smooth_std;
                    
                    fc = fftshift(osl_gauss(st/tres,1,length(dat))');
                    dat2 = fftconv(dat,fc);
                    current_level_results.stdcope(vox,:,c,f,gc)=dat2;
                    
                    %dat2= moving(dat,ceil(2*st/tres));
                    %figure;plot(dat);ho;plot(dat2,'r');
                end;
            end;
            
        end;
    end;
    
    ft_progress('close');
    
    if(current_level.group_varcope_spatial_smooth_fwhm>0 && ~current_level.space_average)
        
        disp(['Do spatial smoothing of group varcopes with FWHM=' num2str(current_level.group_varcope_spatial_smooth_fwhm) ' mm']);
        if strcmp(current_level_results.recon_method,'none'),
            warning('Not implemented for sensor space - no spatial smoothing applied');
        else
                        
            for gc=1:length(group_contrast),
                for f=1:nfreq,
                    dat=permute(current_level_results.stdcope(:,:,c,f,gc),[1,2,3,4,5]);
                    
                    % NOTE - This code is translated as is, but it doesn' t look like it would run - need to check which filename to use for curr_level_mask
                    current_level_results.stdcope(:,:,c,f,gc)=smooth_vol(dat, curr_level_mask_fname, current_level.group_varcope_spatial_smooth_fwhm);

                end;
            end;
                        
        end;
    end;
    
    if(store_lower_level_copes)
        current_level_results.lower_level_copes{c}=cope_smooth_lower_level;
        current_level_results.lower_level_stdcopes{c}=stdcope_smooth_lower_level;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%
%% store results
current_level_results.current_level_contrast_name=current_level.group_contrast_name;
current_level_results.lower_level_contrast_name=lower_level_results.first_level_contrast_name;
current_level_results.current_level=current_level;
current_level_results.lower_level_name=lower_level_results.name;
current_level_results.level=2;
   
%%%%%%%%%%%%%%%%%%%
%% if in sensor space write out stats as SPM MEEG objects
if strcmp(current_level_results.recon_method,'none') && ~current_level.space_average, % sensor space analysis                    

    S4=[];
    S4.oat=oat;       
    S4.stats=current_level_results;
    S4.write_copes=1;
    [current_level_results.D_tstat, current_level_results.D_cope]=oat_save_spm_stats(S4); 

end;

%%%%%%%%%%%%%%%%%%%
%% create report
try
    current_level_results.report = oat_group_level_stats_report(oat,current_level_results);
catch   
end

%%%%%%%%%%%%%%%%%%%
%% save results to file
disp(['Saving statistics in file ' oat.source_recon.dirname '/' current_level_results.fname]);
oat_save_results( oat, current_level_results );

results_fnames=current_level_results.fname;

oat.group_level.results_fnames=results_fnames;

if (strcmp(current_level_results.recon_method,'none') | current_level.space_average),
else
    disp(['To create niftii files from this use a call to oat_save_nii_stats']);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sub funcs
function v2 = smooth_vol(vol_as_matrix, mask_fname, fwhm)

    [mask,res,xform] = nii.load(mask_fname);
    mask(mask>0) = 1; % Turn it into a mask

    sd = fwhm/2.3; % std spatial smoothing, in mm
    sd = sd/res(1); % standard deviation in voxels

    smooth_mask = smooth3(mask,'gaussian',5,sd);
    smooth_mask(~mask) = 0;

    v2 = nan(size(vol_as_matrix));

    for j = 1:size(vol_as_matrix,2)
        v = matrix2vols(vol_as_matrix(:,j),mask);
        v = smooth3(v,'gaussian',5,sd);
        v(~mask) = 0;
        v = v./smooth_mask;
        v2(:,j)=vols2matrix(v,mask);
    end
end

