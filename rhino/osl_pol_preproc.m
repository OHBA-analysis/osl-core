function D=osl_pol_preproc(D,do_plots,remove_nose);
% Function which runs some sanity checks on acquired head points; will
% delete any points which look silly. Optionally also removes points
% (approximately) around the nose.

% If a co-reg has already been run, we give the user the option of doing an
% exploratory second attempt at coregistration excluding any previously
% detected erroneous head points.

% 1) CHECK FOR NONSENSICAL POINTS AND REMOVE THEM
% 2) DELETE POINTS AROUND THE NOSE (OPTIONAL)

% Will one day do:
% 3) DO EXPLORATORY SECOND CO-REG AFTER REMOVING POINTS (IF D OBJECT HAS
%    ALREADY BEEN CO-REGGED) [not yet implemented]

% Inputs (required):
%   D - M/EEG object. If this has been co-registered, some summary
%   statistics will be calculated.

% Inputs (optional):
%   do_plots - boolean. true/false. Whether or not to show logical sphere
%   and nose points being removed (if selected).
%   remove_nose - boolean. true/false. Tries to find points around the nose
%   and removes these.

% Returns:
%   D - M/EEG object, with outliers (and nose points, if selected),
%   removed.

%   orig_hdpnts_for_<D.fname>.mat. A file containing the original head
%   shape points, should you ever need them again. This is stored in the
%   same directory as the D object.

% Example usage:
% D=osl_pol_preproc(D,true,true)

% TO DO:
% Re-run co-reg if user desires, after bad points are removed. User will
% need to specify co-reg options.

% Example data sets:
% D=spm_eeg_load('/Users/rtimms/Documents/DPhil/MEG_UK_Partnership_2019/2020_revisit/sub-astnew001/processed_data/spm_sub01_task-resteyesclosed_meg.mat');
% D=spm_eeg_load('/Users/rtimms/Documents/DPhil/MEG_UK_Partnership_2019/4ryan/sub-oxf001/processed_data/spm_sub01_task-resteyesclosed_meg.mat');

% Ryan Timms (@blobsonthebrain), Andrew Quinn, 2020.

% Check input arguments
%--------------------------------------------------------------------------
try, vis_me = do_plots;   catch, vis_me = true;        end
try, remove_nose;        catch,  remove_nose=false;    end
try, re_run_coreg;        catch, re_run_coreg=false;   end


% Check - has coregistration even been run?
%--------------------------------------------------------------------------
try
    D.inv{1}.mesh.tess_rhino;
    coreg_method='rhino';
catch
    try
        D.inv{1};
        coreg_method='SPM';
    catch
        coreg_method='none';
        warning('Couldn''t detect an existing co-reg. You won''t be able to do any exploratory registrations');
    end
end

% 1) CHECK FOR NONSENSICAL POINTS
% See if there are any points which lie outside of a 15mm radius
% sphere, situated at the origin
%--------------------------------------------------------------------------
pol_points=D.fiducials.pnt;
orig_points=pol_points; % keep a copy of the original head points!
r=150;
fprintf('\n------------------------------------------------\nChecking if all supplied points look sensible\n------------------------------------------------');
outside=(sqrt(sum(pol_points.^2,2))>r);
if sum(outside)==0;
    fprintf('\nAll head points look okay!\n');
    points_to_del=[];
else
    points_to_del=find(outside);
    fprintf('\nDetected points that are unlikely to be on the scalp. Deleting the following points: %0.0f\n',points_to_del)
    pol_points(points_to_del,:)=[];
end
if vis_me;
    figure;
    hold all;
    scatter3(pol_points(:,1),pol_points(:,2),pol_points(:,3),'bo','filled','MarkerEdgeColor','k')
    scatter3(pol_points(points_to_del,1),pol_points(points_to_del,2),pol_points(points_to_del,3),'ro','filled','MarkerEdgeColor','k')
    axis equal
    axis off
    legend('Okay points','Point(s) which will be removed');
    [X,Y,Z] = sphere;
    surf(r*X,r*Y,r*Z,'FaceAlpha',0.2)
    colormap gray
    drawnow
end
    
% Delete the trashy head points
R=D.fiducials;
R.label(points_to_del)=[];
R.pnt(points_to_del,:)=[];
D=fiducials(D,R);
warning('DELETING EXISTING HEAD POINTS!');
warning('THE ORIGINAL POINTS ARE SAVED AT THE FOLLOWING LOCATION:');
orig_points_fname=[D.path,'/orig_hdpnts_for_',D.fname];
fprintf(orig_points_fname)
fprintf('\n');
save(orig_points_fname,'orig_points');
D.save;


% 2) DELETE POINTS AROUND THE NOSE (OPTIONAL)
if remove_nose;
    fprintf('\n------------------------------------------------\nRemoving points from nose\n------------------------------------------------\n');
    
    pnt = D.fiducials.pnt;
    labz = D.fiducials.fid.label;
    % Find what index corresponds to the nose...
    nose_ind=find(ismember(labz,{'Nasion','NASION','Nas','nas'}));
    
    nas = D.fiducials.fid.pnt(nose_ind,:);
    nas(3) = nas(3) - 40; % drop nasion by ~4cm
    distances = sqrt((nas(1)-pnt(:,1)).^2 + (nas(2)-pnt(:,2)).^2 + (nas(3)-pnt(:,3)).^2 );
    keeps = distances>60; % remove points within 6cm
    fids = D.fiducials;
    fids.pnt = fids.pnt(keeps,:);
    if vis_me;
        figure;
        hold all;
        scatter3(pnt(:,1),pnt(:,2),pnt(:,3),'ro','filled','MarkerEdgeColor','k')
        scatter3(fids.pnt(:,1),fids.pnt(:,2),fids.pnt(:,3),'bo','filled','MarkerEdgeColor','k')
        axis image
        view([-113 30])
        legend('Removed points','Maintained points');
    end
    fprintf('Done.\n');
    D = fiducials(D,fids);
    D.save % Update D object's head points
end


% 3) DO EXPLORATORY SECOND CO-REG AFTER REMOVING POINTS (IF D OBJECT HAS
%    ALREADY BEEN CO-REGGED)
if strcmp(coreg_method,'none')~=1;
    fprintf('\n------------------------------------------------\nCalculating polhemus error statistics\n------------------------------------------------\n');
    
    % Calculate existing polhemus error
    %--------------------------------------------------------------------------
    if strcmp(coreg_method,'rhino')==1;
        mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);
        err = rhino_cperror(mesh_rhino.vert,pol_points);
    elseif strcmp(coreg_method,'SPM')==1;
        mesh = spm_eeg_inv_transform_mesh(D.inv{1}.datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
        scalp=mesh.tess_scalp;
        err = rhino_cperror(scalp.vert,pol_points);
    end
    
    orig_mean_err=mean(err);
    
    fprintf('Running outlier detection...');
    if find(isoutlier(err))>0;
        fprintf('\nDetected outlying head points.\n');
        fprintf('We''ll temporarily remove these points, and then re-run the coreg to see if this improves things');
        
        points_to_del=[points_to_del; find(isoutlier(err))];
    else
        fprintf('\nCouldn''t detect any outlying head points.\n');
    end
    
    % Re-Calculate polhemus error, without the crap pol points
    pol_points(points_to_del,:)=[];
    if strcmp(coreg_method,'rhino')==1;
        mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);
        err = rhino_cperror(mesh_rhino.vert,pol_points);
    elseif strcmp(coreg_method,'SPM')==1;
        scalp=mesh.tess_scalp;
        err = rhino_cperror(scalp.vert,pol_points);
    else
        error('Unbeknown error');
    end
    
    orig_mean_err_excluding_outliers=mean(err);
    fprintf('\n\nInitial mean error - all points: %0.2fmm\n',mean(orig_mean_err))
    fprintf('Initial mean error - excluding detected outliers: %0.2fmm\n',mean(orig_mean_err_excluding_outliers))
end

D=spm_eeg_load(D.fullfile);

% THE FOLLOWING HAS NOT BEEN TESTED AND WILL BE INCLUDED IN A FUTURE
% ITERATION OF THE CODE

% if re_run_coreg; % Re-run coreg
%     fprintf('\n------------------------------------------------\nRe-running co-reg!\n------------------------------------------------\n');
%     
%     % Make a copy of the original D object
%     copy_of_D=copy(D,'copy_of_D');
%     
%     % Delete the bad head points
%     R=copy_of_D.fiducials;
%     R.label(points_to_del)=[];
%     R.pnt(points_to_del,:)=[];
%     copy_of_D=fiducials(copy_of_D,R);
%     copy_of_D.save;
%     copy_of_D=spm_eeg_load(copy_of_D.fullfile);
%     
%     
%     coreg_settings = struct;
%     coreg_settings.D = copy_of_D.fullfile;
%     coreg_settings.mri =copy_of_D.inv{1}.mesh.sMRI;
%     coreg_settings.useheadshape = true;
%     coreg_settings.forward_meg = 'Single Sphere';
%     coreg_settings.use_rhino = true;
%     coreg_settings.fid.label.nasion='Nasion';
%     coreg_settings.fid.label.lpa='LPA';
%     coreg_settings.fid.label.rpa='RPA';
%     copy_of_D = osl_headmodel(coreg_settings);
%     
%     
%     % Re-compute stats
%     pol_points=copy_of_D.fiducials.pnt;
%     
%     if strcmp(coreg_method,'none')~=1;
%         if strcmp(coreg_method,'rhino')==1;
%             mesh_rhino = spm_eeg_inv_transform_mesh(copy_of_D.inv{1}.datareg(1).fromMNI*copy_of_D.inv{1}.mesh.Affine, copy_of_D.inv{1}.mesh.tess_rhino);
%             err = rhino_cperror(mesh_rhino.vert,pol_points);
%         elseif strcmp(coreg_method,'SPM')==1;
%             scalp=mesh.tess_scalp;
%             err = rhino_cperror(scalp.vert,pol_points);
%         end
%     end
%     fprintf('\nFinal mean error, after second co-reg: %0.2fmm\n',mean(err))
% end

end