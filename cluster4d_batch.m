function cluster_stats=cluster4d_batch(S)

% cluster_stats=cluster4d_batch(S)
%
% batch script to run 4D-permutation test on MEG data
% requires input images for each timepoint stored in a single directory:%
%
% - INPUT arguments:
% S.dirname is folder containing a 4D (voxels*subjects) image called
%       'all_subsXXXX.nii.gz', where XXXX is each timepoint in S.tp
%  .tp is an array of all timepoints
%  .nP is number of permutations (must be multiple of 100, default 5000)
%  .write_cluster_script, just writes out generated scripts for running on a
%  UNIX cluster, otherwise will run a single thread on the host comptuer
%       
%  .thresh is T-statistic threshold to apply to images 
%
% - OUTPUT directories are in S.dirname/permutations
%
% Laurence Hunt and Tom Nichols and Mark Woolrich, February 2010/11 (beta
% version!)

try, dirname = S.dirname; catch, error('Must specify S.dirname'); end
try, tp = S.tp; catch, error('Must specify S.tp'); end
try, times = S.times; catch, error('Must specify S.times'); end
try, nP = S.nP; catch, nP = 5000; end %number of permutations
try, write_cluster_script = S.write_cluster_script; catch, write_cluster_script = 1; end
try, thresh = S.thresh; catch, thresh = 2.3; end
try, mask_fname = S.mask_fname; catch, mask_fname='/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'; end
try, gridstep = S.gridstep; catch, error('Must specify S.gridstep'); end
try, contrasts = S.contrasts; catch, error('Must specify S.contrasts'); end
try, X = S.X; catch, error('Must specify S.X'); end
try, subdirname=S.subdirname; catch, subdirname='permutations'; end;
try, S.group_varcope_spatial_smooth_std=S.group_varcope_spatial_smooth_std; catch S.group_varcope_spatial_smooth_std=0; end;
try, fsl_version_4p1=S.fsl_version_4p1; catch fsl_version_4p1=1; end;
try, S.matlab_exe_name=S.matlab_exe_name; catch S.matlab_exe_name='matlab'; end;

if(rem(nP,100)~=0)
    error('S.nP must be divisible by 100.');
end


cmd = ['mkdir ' dirname '/' subdirname];
unix(cmd);
scriptdir = [dirname '/' subdirname '/clusterscripts'];
mkdir(scriptdir);
lfname = [scriptdir '/logfiles'];
mkdir(lfname);

delete([scriptdir '/DO.sh']);

fidM =fopen([scriptdir '/DO.sh'],'w'); %overarching script
fprintf(fidM,'jid3=1\n');

% save design matrix and contrasts
save_vest(X',[scriptdir '/design.mat']);
save_vest((contrasts)',[scriptdir '/design.con']);    

for r = 1:nP/100; 
    %% 5.1 -run randomise
    fid = fopen([scriptdir '/split_merge_randomise5_' num2str(r) '_1.sh'],'w');
    for t = tp;
        permdir = sprintf('%s/%s/%04.0f',dirname,subdirname,t);
        mkdir(permdir);
        fprintf(fid,['randomise -d ' [scriptdir '/design.mat'] ' -t ' [scriptdir '/design.con'] ' -i %s/allsubs_time%04.0f -o %s/stats%04.0f -c %1.2f -R -n 100 ' ...
            '--seed=%0.0f -v ' num2str(S.group_varcope_spatial_smooth_std) ' --permout -m %s/mask_time%04.0f \n'],...
            dirname,t,permdir,t,thresh,r,dirname,t); % -c means cluster-based thresholding
    end
    fclose(fid);
    
    %% 5.2 - produce clusters
    fid = fopen([scriptdir '/split_merge_randomise5_' num2str(r) '_2.sh'],'w');
    for t = tp;
        permdir = sprintf('%s/%s/%04.0f',dirname,subdirname,t);
        for p = 1:100;
            if fsl_version_4p1
                fprintf(fid,['cluster -i %s/stats%04.0f_rawstat_tstat1_%05.0f -t %1.2f -o %s/cindex%05.0f; '],...
                    permdir,t,p,thresh,permdir,(r-1)*100+p); % -c means cluster-based thresholding
            else
                fprintf(fid,['cluster -i %s/stats%04.0f_vox_tstat1_perm%05.0f -t %1.2f -o %s/cindex%05.0f; '],...
                permdir,t,p,thresh,permdir,(r-1)*100+p); % -c means cluster-based thresholding
            end;
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    %% 5.3 - delete randomise output
    fid = fopen([scriptdir '/split_merge_randomise5_' num2str(r) '_3.sh'],'w');
    for t = tp;
        permdir = sprintf('%s/%s/%04.0f',dirname,subdirname,t);
        fprintf(fid,'rm ');
        for p = 1:100
            if fsl_version_4p1
                fprintf(fid,'%s/stats%04.0f_rawstat_tstat1_%05.0f.nii.gz ',permdir,t,p);
            else
                fprintf(fid,'%s/stats%04.0f_*_tstat1_perm%05.0f.nii.gz ',permdir,t,p);
            end;
        
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    %% write into overarching script
    if(~write_cluster_script),
        fprintf(fidM,'%s/split_merge_randomise5_%0.0f_1.sh  \n',scriptdir,r);
        fprintf(fidM,'%s/split_merge_randomise5_%0.0f_2.sh  \n',scriptdir,r);
        fprintf(fidM,'%s/split_merge_randomise5_%0.0f_3.sh  \n',scriptdir,r);
    else
        fprintf(fidM,'jid1=`fsl_sub -q veryshort.q -j $jid3 -t %s/split_merge_randomise5_%0.0f_1.sh -l %s` \n',scriptdir,r,lfname);
        fprintf(fidM,'jid2=`fsl_sub -q veryshort.q -j $jid1 -t %s/split_merge_randomise5_%0.0f_2.sh -l %s` \n',scriptdir,r,lfname);
        fprintf(fidM,'jid3=`fsl_sub -q veryshort.q -j $jid2 -t %s/split_merge_randomise5_%0.0f_3.sh -l %s` \n',scriptdir,r,lfname);
    end;
end

%% write a matlab call to cluster4d_dist

fsl_sub_rsize=100;

for r = 1:nP/fsl_sub_rsize;
    distfile_save=sprintf('%s/%s/dist%0.0f.mat',dirname,subdirname,r);
    mfile = sprintf('%s/c4ddist_%0.0f.m',scriptdir,r);
    fid = fopen(mfile,'w');
    
    tmp=which('cluster4d_dist');
    if strcmp(tmp,''), error('cluster4d_dist not in path'); end;
    
    fprintf(fid,['addpath(''%s'');\n addpath(''%s/etc/matlab'');\n' ...
        'S.dirname = ''%s''; S.subdirname = ''%s''; S.distfile_save = ''%s''; S.np = %0.0f:%0.0f;'],...
        tmp(1:end-17),getenv('FSLDIR'),dirname,subdirname,distfile_save,(r-1)*fsl_sub_rsize+1,r*fsl_sub_rsize);
    fprintf(fid,['S.tp = ',mat2str(tp) ';']);
    fprintf(fid,['S.times = ',mat2str(times) ';']);
    fprintf(fid,'S.save_images = 0;\n S.gridstep = %0.0f; \n cluster4d_dist(S);',gridstep);
    fclose(fid);

    if(~write_cluster_script),
        fprintf(fidM,'%s -nojvm -nodisplay -nosplash -singleCompThread \\< %s \n',S.matlab_exe_name,mfile);
    else,
        fprintf(fidM,'jid4=`fsl_sub -q long.q -j $jid3 -l %s %s -nojvm -nodisplay -nosplash -singleCompThread \\< %s > %s/log_cluster4d_dist.log  2>&1` \n',lfname,S.matlab_exe_name,mfile,lfname);
    end;
      
end


%% write a matlab call to combine distributions together

mfile = sprintf('%s/combine_dists.m',scriptdir);
fid = fopen(mfile,'w');

fprintf(fid,['clusterstats_tmp.dist=[];\n']);

for r = 1:nP/fsl_sub_rsize;

    distfile_save=sprintf('%s/%s/dist%0.0f.mat',dirname,subdirname,r);

    fprintf(fid,['clusterstats=load(''' distfile_save ''');\n']);
    fprintf(fid,['clusterstats_tmp.dist=[clusterstats_tmp.dist; clusterstats.dist];\n']);
    
end;

fprintf(fid,['clusterstats.dist=clusterstats_tmp.dist;']);
distfile_save=sprintf('%s/%s/dist_combined.mat',dirname,subdirname);
fprintf(fid,['save(''' distfile_save ''', ''-struct'', ''clusterstats'');']);
      
fclose(fid);
    
if(~write_cluster_script),
    fprintf(fidM,'%s -nojvm -nodisplay -nosplash -singleCompThread \\< %s \n',S.matlab_exe_name,mfile);
else,
    fprintf(fidM,'fsl_sub -q long.q -j $jid4 -l %s %s -nojvm -nodisplay -nosplash -singleCompThread \\< %s > %s/log_cluster4d_dist.log  2>&1 \n',lfname,S.matlab_exe_name,mfile,lfname);
end;
    
fclose(fidM);
cmd = ['chmod -R u+x ' scriptdir];
unix(cmd);
cmd = [scriptdir '/DO.sh'];

%% run, if on cluster
cluster_stats=[];

if(~write_cluster_script),
            
    disp(['Script is autorunning on host system:']);
    disp(cmd);
    s = unix(cmd);

    disp(['Script finished.']);    
    disp('Script return code:');
    disp(s);    
    
    try,
        cluster_stats=load(distfile_save);
    catch,
        disp(['Cluster 4D script has failed. Look in the log files in dir ' lfname]);
    end;
    
else,
    disp(['Command to run is: ' cmd]);
    disp(['Logfiles in: ' scriptdir '/logfiles']);
    disp(['When finished, cluster_stats results will be in: cluster_stats=load(''' distfile_save ''');']);
end

