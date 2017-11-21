function res = osl_render_vols_to_surf(S)

% res = osl_render_vols_to_surf(S)
%
% For S.type='freesurfer' you need Freesurfer and ImageMagick installed
% For S.type='workbench' you need HCP workbench installed
%
% E.g.:
% S=[];
% S.time_indices=[1];S.vol='/Users/woolrich/homedir/matlab/osl_testdata_dir/ctf_fingertap_subject1_data/subj1_results_beta.oat/wholebrainnew_stats_dir/tstat3_2mm';
% S.time_indices=[60:70];S.vol='/Users/woolrich/homedir/matlab/osl_testdata_dir/faces_subject1_data/spm8_meg1.mat_wideband_beamform.oat/session1_wholebrain_first_level_dir/tstat3_2mm.nii.gz'
% S.view=1;  S.minmax=[3, 7]; S.type='freesurfer'; 
%
% Mark Woolrich 2013

OSLDIR = getenv('OSLDIR');

try, time_indices=S.time_indices; catch S.time_indices=1; end;

try, time_size=S.time_indices(end)-S.time_indices(1)+1; catch time_size=1; end;
try, time_lower=S.time_indices(1)-1; catch time_lower=0; end;

switch S.type,
    case {'freesurfer'},
    
        if(isempty(getenv('FREESURFER'))),
            warning('Need freesurfer installed and with environmental variables set');         

            FREESURFER_HOME='/Applications/freesurfer';
            SUBJECTS_DIR=[OSLDIR '/std_masks/fs_subjects'];
            setenv('SUBJECTS_DIR', SUBJECTS_DIR);
            setenv('FREESURFER_HOME', FREESURFER_HOME);
        else
            FREESURFER_HOME=getenv('FREESURFER');
            SUBJECTS_DIR=getenv('SUBJECTS_DIR');
        end;        

        workingdir=[S.vol '_dir'];
        mkdir(workingdir);    
        cd(workingdir);
        
        runcmd(['rm -rf ' workingdir '/volumes']);
        mkdir([workingdir '/volumes']);     
        
        runcmd(['rm -rf ' workingdir '/surface_frames']);
        mkdir([workingdir '/surface_frames']);
            
        times=time_lower:(time_lower+time_size-1);

        for ind=1:length(times),
            ind/length(times)
            time=times(ind);

            filename=[workingdir '/volumes/tfile_vol' num2str(time) '.nii.gz'];            

            runcmd(['fslroi ' S.vol ' ' filename ' ' num2str(time) ' ' num2str(1)]);
        
            % for left hemi:
            runcmd([FREESURFER_HOME '/bin/mri_vol2surf --src ' filename ' --src_type nifti --srcreg ' OSLDIR '/std_masks/MNI_to_struct.dat --hemi lh --projfrac 0.5 --out ' workingdir '/surface_frames/TS' num2str(time) '_vol-lh.w --out_type paint']);
            runcmd([FREESURFER_HOME '/bin/tksurfer FS_spm_CanonicalBrain lh inflated -curv -gray -overlay surface_frames/TS' num2str(time) '_vol-lh.w -overlay-reg-identity -fminmax ' num2str(S.minmax(1)) ' ' num2str(S.minmax(2)) ' -colscalebarflag 1 -tcl ' OSLDIR '/osl/render_tcl.tcl']);
            runcmd(['cp ' workingdir '/surface_frames/TS' num2str(time) '_vol-lh.w_all.tiff ' workingdir '/volumes']);

            % for right hemi:
            runcmd([FREESURFER_HOME '/bin/mri_vol2surf --src ' filename ' --src_type nifti --srcreg ' OSLDIR '/std_masks/MNI_to_struct.dat --hemi rh --projfrac 0.5 --out ' workingdir '/surface_frames/TS' num2str(time) '_vol-rh.w --out_type paint']);
            runcmd([FREESURFER_HOME '/bin/tksurfer FS_spm_CanonicalBrain rh inflated -curv -gray -overlay surface_frames/TS' num2str(time) '_vol-rh.w -overlay-reg-identity -fminmax ' num2str(S.minmax(1)) ' ' num2str(S.minmax(2)) ' -colscalebarflag 1 -tcl ' OSLDIR '/osl/render_tcl.tcl']);                    
            runcmd(['cp ' workingdir '/surface_frames/TS' num2str(time) '_vol-rh.w_all.tiff ' workingdir '/volumes']);

            res.left_tiff_names{ind}=[workingdir '/volumes/TS' num2str(time) '_vol-lh.w_all.tiff'];        
            res.right_tiff_names{ind}=[workingdir '/volumes/TS' num2str(time) '_vol-rh.w_all.tiff'];
        end;
        
        %left
        cmd=['open ' res.left_tiff_names{ind} ' &'];
            
        disp(cmd);
        
        if(S.view)
            runcmd(cmd);
        end;
        
        %right
        cmd=['open ' res.right_tiff_names{ind} ' &'];
            
        disp(cmd);
        
        if(S.view)
            runcmd(cmd);
        end;                
    
        res.workingdir=workingdir;
        res.times=times;

    case {'workbench'},        
        
        if(time_size>1)
            error('Workbench only works with single volumes - try the freesurfer option');
        end;
               
        workbench_binaries_dir='/Users/woolrich/homedir/hcp/workbench/bin_macosx64';

        output_left=[S.vol '_left.func'];
        output_right=[S.vol '_right.func'];

        surf_left=[OSLDIR '/std_masks/ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii'];
        surf_right=[OSLDIR '/std_masks/ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii'];

        runcmd(['surf_proj --data=' S.vol ' --surf=' surf_left ' --meshref=' S.vol ' --out=' output_left ' --surfout ']);
        runcmd(['surf_proj --data=' S.vol ' --surf=' surf_right ' --meshref=' S.vol ' --out=' output_right ' --surfout ']);

        output_right=[output_right '.gii'];

        output_left=[output_left '.gii'];

        runcmd([workbench_binaries_dir '/wb_command -set-structure ' output_right ' CORTEX_RIGHT']);
        runcmd([workbench_binaries_dir '/wb_command -set-structure ' output_left ' CORTEX_LEFT']);

        cmd=[workbench_binaries_dir '/workbench ' surf_left ' ' surf_right ' ' output_left ' ' output_right ' &'];
        disp(cmd);
        
        if(S.view)
            runcmd(cmd);
        end;
       
end;

