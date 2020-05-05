function osl_RAC(D,varargin);
% OSL's Rhino's Automatic Co-registration Check.

% Function which takes in a co-regged D object and runs some simple sanity
% checks to assist the user with assessing whether or not a co-registration
% has been successful, or at least plausibly so. Makes use of the inhull
% and findTriMeshHoles functions from the MATLAB file exchange. A short and
% long form report are output as .txt files to the output_dir. A co-reg
% visualisation Figure is also saved to the output_dir, along with a
% histogram of digitised head point to scalp distances. A "traffic_report"
% traffic-light style Figure is also produced showing whether or not a
% sanity check has been passed (green), plausibly passed (orange) or failed
% (red). A white traffic light means that we were not able to run a check
% (e.g. in the absence of digitised head points).

% This function can be particularly useful when dealing with datasets which
% contain a large number of subjects.

% Should now support M/EEG fusion D objects, but note that only the MEG
% co-reg is reported. Also returns STD of co-reg digitised head point
% error. Now reporting both the RMS and dx/y/z error to check for
% systematic errors

% Example usage:
% osl_RAC(D,'output_dir') % to save in 'output_dir'
% osl_RAC(D) % to save report in pwd

% Ryan Timms, 2020.

% References:
%
% John D'Errico (2020). Inhull
% (https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull),
% MATLAB Central File Exchange. Retrieved May 5, 2020.
%
% Audrey Cheong (2020). Find Holes in Triangular Mesh
% (https://www.mathworks.com/matlabcentral/fileexchange/62419-find-holes-in-triangular-mesh),
% MATLAB Central File Exchange. Retrieved May 5, 2020.


% An example CTF dataset
% D=spm_eeg_load('/Users/rtimms/Documents/DPhil/Source_recon/digit_tap/pre_processed_data/08966/fedff08966080_Ellie_20160210_01.mat');

% Example Elekta datasets
% D=spm_eeg_load('/Users/rtimms/Documents/DPhil/MEG_UK_Partnership_2019/2020_revisit/sub-astnew001/processed_data/spm_sub01_task-resteyesclosed_meg.mat');
% D=spm_eeg_load('/Users/rtimms/Documents/DPhil/MEG_UK_Partnership_2019/4ryan/sub-oxf001/processed_data/spm_sub01_task-resteyesclosed_meg.mat');

fprintf('==================================================================')
fprintf('\nCalling out the RAC...\n');
fprintf('==================================================================\n')

if nargin<2;
    output_dir=string(pwd);
else
    output_dir=varargin(1);
    mkdir(string(output_dir))
end

% Has a co-reg been run?
try
    D.inv{1}.mesh.tess_rhino;
    coreg_method='rhino';
catch
    try
        D.inv{1};
        coreg_method='SPM';
    catch
        error('Neither RHINO or SPM co-reg has been run! Please check the D object');
    end
end

% Need to check if the data set has EEG channels in it, too. Run a quick
% check now before doing anything else.
if strcmp(D.inv{1}.datareg(1).modality,'EEG')==1;
    if strcmp(D.inv{1}.datareg(2).modality,'MEG')==1;
        MEG_mod_ind=2;
    end
elseif strcmp(D.inv{1}.datareg(1).modality,'MEG')==1;
    MEG_mod_ind=1;
else
    error('Couldn''t detect any MEG channels in D object. EEG co-reg check not supported');
end



if strcmp(D.modality,'Multimodal')==1; % Thanks, CBU
    scanner_type=D.inv{1}.datareg(MEG_mod_ind).sensors.coordsys;
    multimodal=1;
else
    % Elekta or CTF data?
    scanner_type=D.inv{1}.datareg(MEG_mod_ind).sensors.coordsys;
    multimodal=0;
end
if multimodal==1;
    warning('Detected a fused M/EEG dataset. Are you from Cambridge? The RAC hasn''t been tested too much with these kind of data. Proceed with caution.');
end

% Check 1: did all of the fiducial points end up inside the dewar? n.b.
% that we could have some polhemus points outside of the dewar (e.g. on the
% nose/cheek), so best just to raise a warning if we see *loads* of these
% outside of the head.

if strcmp(scanner_type,'ctf')==1;
    % For CTF data, get all of the gradiometer positions
    gradiometer_indices=find(strcmp(D.chantype,{'MEGGRAD'}));
elseif strcmp(scanner_type,'neuromag')==1;
    if multimodal==0;
        % For Elekta data, get all of the MEG channel locations
        gradiometer_indices=[1:length(D.inv{1}.datareg.sensors.chanpos)];
    else
        gradiometer_indices=[1:length(D.inv{1}.datareg(MEG_mod_ind).sensors.chanpos)];
    end
end

if multimodal==0;
    gradiometer_locations=D.inv{1}.datareg.sensors.chanpos(gradiometer_indices,:);
else
    gradiometer_locations=D.inv{1}.datareg(MEG_mod_ind).sensors.chanpos(gradiometer_indices,:);
end

fprintf('\n==================================================================\n')
fprintf('CHECK 1: Did the fiducials end up inside the dewar?\n\n')

% Now find the extremeties of the dewar; max/min x/y/z location
[max_x,min_x,max_y,min_y,max_z,min_z]=dewar_extremities(gradiometer_locations);

% Run fiducial sanity check
if multimodal==0;
    MEG_fids=D.inv{1}.datareg.fid_eeg.fid;
else
    MEG_fids=D.inv{1}.datareg(MEG_mod_ind).fid_eeg.fid;
end
fid_in_score=0;
str_log='';
for i=1:length(MEG_fids.label);
    fprintf('\nChecking that the "%s" is inside the dewar...',MEG_fids.label{i})
    current_fid=MEG_fids.pnt(i,:);
    check_x=current_fid(1)<max_x && current_fid(1)>min_x;
    check_y=current_fid(2)<max_y && current_fid(2)>min_y;
    check_z=current_fid(3)<max_z && current_fid(3)>min_z;
    if check_x+check_y+check_z==3;
        fprintf(' it looks like it!\n');
        fid_in_score=fid_in_score+1;
    else
        fprintf(' it doesn''t appear to be. Will take a note of this.\n\n');
        str_log=strcat(str_log,MEG_fids.label{i},{', '})
        if ismember(MEG_fids.label{i},{'Nasion','Nas'})==1;
            fprintf('(Note that for the nasion this could be a feasible co-reg result.\n)');
        end
    end
end
if fid_in_score==3;
    fprintf('\nGood news - looks like all fiducial points are within the dewar\n');
elseif fid_in_score<3 && fid_in_score>0;
    fprintf('Hmm. Looks like at least one of the fiducial markers ended up outside of the scanner.\nWe recommend that you take a look again at your co-reg\n');
elseif fid_in_score==0;
    fprintf('None of the fiducials ended up inside the dewar. Something bad may have happened.\n');
    fid_in_score=0; %set fiducial score to be equal to 0 - we'll come back to this at the end
end

% Tidy up
clear check_x check_y check_z current_fid max_x min_x max_y min_y max_z min_z where who

% Check 2: the L/RPA points should be roughly on the same plane between the
% ears as one another, i.e. on a line through the z axis. Check that this
% is the case! Note that at least for the SPM co-reg, the M/EEG fiducials
% define the origin, hence they will necessarily be on the same z. This is
% why we use the datareg.fid_mri.fid.pnt points

fprintf('\n==================================================================\n')
fprintf('CHECK 2: Are the L/RPA points aligned after co-reg?\n\n')

% First we need to find out what axis is used to define the front of the
% head. Easily achieved by finding the maximum value in the nasion fiducial
% position. Have to account for different naming conventions across
% sites/scanners
tmp_nas=find(ismember(MEG_fids.label,{'Nasion','NASION','Nas','nas'}));
if sum(tmp_nas)==0;
    error('Couldn''t detect a nasion point, sorry.');
end
[where,who]=max(abs(MEG_fids.pnt(tmp_nas,:)));

store_points=[];
% Extract just the L/RPA points:

if multimodal==0;
    for i=1:length(D.inv{1}.datareg.fid_mri.fid.label);
        if strcmpi(D.inv{1}.datareg.fid_mri.fid.label{i},'rpa') || strcmpi(D.inv{1}.datareg.fid_mri.fid.label{i},'lpa')
            store_points=[store_points; D.inv{1}.datareg.fid_mri.fid.pnt(i,:)];
        end
    end
else %multimodal case. Assume that the MEG fiducials are in datareg(MEG_mod_ind)
    for i=1:length(D.inv{1}.datareg(MEG_mod_ind).fid_mri.fid.label);
        if strcmpi(D.inv{1}.datareg(MEG_mod_ind).fid_mri.fid.label{i},'rpa') || strcmpi(D.inv{1}.datareg(MEG_mod_ind).fid_mri.fid.label{i},'lpa')
            store_points=[store_points; D.inv{1}.datareg(MEG_mod_ind).fid_mri.fid.pnt(i,:)];
        end
    end
    
end
vert_misalign=abs(store_points(1,who)-store_points(2,who));
clear store_points
fprintf('Vertical misalignment is %f mm;',vert_misalign)
if vert_misalign<10;
    fprintf(' that''s not too bad.\n');
    LR_vert_score=1;
else
    fprintf(' that looks a bit suspect. Will take a note.\n');
    LR_vert_score=0;
    if exist('str_log','var') == 1;
        str_log=strcat(str_log,'. L/RPA points didn''t end up on the same plane');
    else
        str_log='L/RPA points didn''t end up on the same plane';
    end
end

fprintf('\n==================================================================\n')
fprintf('CHECK 3: Did the fiducials end up within ±5mm of one another?\n\n')
if multimodal==0;
    MR_fids= D.inv{1}.datareg.fid_mri.fid;
else
    MR_fids= D.inv{1}.datareg(MEG_mod_ind).fid_mri.fid;
end
fid_delta_score=0;
delta_store=[];
fid_label_store={};

% This part needs to have some intelligent way of checking that we're
% comparing nas with nas, lpa w/ lpa and rpa w/ rpa. Unfortunately,
% different centres have different naming conventions. Need to think about
% this. See commented out lines
for i=1:3;
    %     if strcmp(MR_fids.label{i},MEG_fids.label{i})==1;
    delta_store(i,:)=norm(MR_fids.pnt(i,:)-MEG_fids.pnt(i,:));
    diff(i,:)=MR_fids.pnt(i,:)-MEG_fids.pnt(i,:);
    fid_label_store{i}=MR_fids.label{i};
    if norm(MR_fids.pnt(i,:)-MEG_fids.pnt(i,:))<5;
        fid_delta_score=fid_delta_score+1;
    end
    %     else
    %         error('Not comparing like-with-like fiducials. Eek.')
    %     end
end
if fid_delta_score==3;
    fprintf('\nThe co-reg error between MEG & MRI fiducials was sufficiently small.\nRMS error (mm):\n%s: %f  %s: %f  %s: %f',string(fid_label_store(1)),delta_store(1),string(fid_label_store(2)),delta_store(2),string(fid_label_store(3)),delta_store(3));
else
    fprintf('\nHigh co-reg error between MEG & MRI fiducials. Will take a note of this.\nRMS error (mm):\n\n%s: %f  %s: %f  %s: %f',string(fid_label_store(1)),delta_store(1),string(fid_label_store(2)),delta_store(2),string(fid_label_store(3)),delta_store(3));
end
fprintf('\n\nElement wise error (delta x, delta y, delta z):\n');
for i=1:numel(fid_label_store);
    fprintf('\n---------------\n%s:\n---------------\n%2.2fmm %2.2fmm %2.2fmm',string(fid_label_store(i)),diff(i,1),diff(i,2),diff(i,3));
end

fprintf('\n==================================================================\n')
fprintf('CHECK 4: Is the brain inside the skull?\n\n')
if multimodal==0;
    mesh = spm_eeg_inv_transform_mesh(D.inv{1}.datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
else
    mesh = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(MEG_mod_ind).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
end
Mcortex = mesh.tess_ctx;
Miskull = mesh.tess_iskull;

cort_vert    = Mcortex.vert; % cortical mesh
iskull_vert    = Miskull.vert; % inner skull

brain_in_skull_score=0;
check_is_in=inhull(cort_vert,iskull_vert);

% figure;
% scatter3(cort_vert(:,1),cort_vert(:,2),cort_vert(:,3));
% hold all
% scatter3(iskull_vert(:,1),iskull_vert(:,2),iskull_vert(:,3));

% When using custom meshes (e.g. from FreeSurfer output), the cortex can
% sometimes clip in the inner skull. This in practice isn't a huge problem
% (according to Jimmy). Should raise a warning rather than a hard fail
% (i.e. orange rather than red) if this is the case.

% TO DO: automatically check if a custom inner mesh has been used and
% adjust output accordingly
if sum(check_is_in)==length(cort_vert)
    fprintf('\nYes - all of the brain is inside the skull. Phew!\n');
    brain_in_skull_score=2;

elseif sum(check_is_in)>round(0.97*length(cort_vert));
    fprintf('\nAt least 97 percent of the brain was found inside the skull, but some was outside.\nThis can sometimes happen if you used a custom mesh.\n');
    brain_in_skull_score=1;
else
    fprintf('\nUh-oh. The brain is outside of the skull. Better call Mark Woolrich.\n')
end

fprintf('\n==================================================================\n')
fprintf('CHECK 5: Is the skull inside the dewar?\n\n')
% n.b., we can allow some of the skull to be outside of the dewar in the
% vertical direction (i.e. bottom of the skull can be below the dewar).
% This means that we can't use this diagnostic exclusively to say whether
% or not the co-reg worked/failed, but it can definitely be used as a red
% flag!
check_is_in=inhull(iskull_vert,gradiometer_locations);

if sum(check_is_in)==length(iskull_vert);
    fprintf('It looks like all of the skull is inside the dewar.\n');
    skull_in_dewar_score=2;
elseif sum(check_is_in)>=round(0.80*length(iskull_vert));
    fprintf('The majority of the skull is inside the dewar, but some points are outside it.\n');
    skull_in_dewar_score=1;
else
    fprintf('A significant number of skull points are outside of the dewar. Something could be seriously wrong.\n');
    skull_in_dewar_score=0;
end

fprintf('\n==================================================================\n')
fprintf('CHECK 6: Were any head points used in the co-reg?\n\n')
if multimodal==0
    pol_points=D.inv{1}.datareg.fid_eeg.pnt;
else
    pol_points=D.inv{1}.datareg(MEG_mod_ind).fid_eeg.pnt;
end
pol_point_score=[];
if numel(pol_points)>9 % could be [LPA, RPA, Nas x3]
    fprintf('Yes. Good!\n');
    pol_point_score=1;
else
    fprintf('No digitised head points detected.\n') % could be headcast data?
    pol_point_score=0;
end

fprintf('\n==================================================================\n')
fprintf('CHECK 7: Was a custom MRI used to co-reg?\n\n')
custom_anat_score=[];
if D.inv{1}.mesh.template==0;
    fprintf('Yes. Good!\n');
    custom_anat_score=1;
else
    fprintf('No. Taking a note of this.\n');
    custom_anat_score=0;
end

fprintf('\n==================================================================\n')
fprintf('CHECK 8: Polhemus point statistics\n\n')
if pol_point_score==0;
    fprintf('Cannot do any head point statistics without any head points!\n');
    fprintf('Skpping this check.\n');
    err=NaN; %set to NaN
else
    if strcmp(coreg_method,'rhino')==1;
        mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);
        err = rhino_cperror(mesh_rhino.vert,pol_points);
    elseif strcmp(coreg_method,'SPM')==1;
        scalp=mesh.tess_scalp;
        err = rhino_cperror(scalp.vert,pol_points);
    else
        error('Unbeknown error');
    end
    
    n_hist_pts=100;
    [n,x]=hist(err,n_hist_pts);
    
    
    % Create green to yellow to red colormap
    G2R=[linspace(0.5,1,n_hist_pts);linspace(1,0,n_hist_pts);zeros(1,n_hist_pts)];
    G2R=G2R';
    
    figure('Position',[10 10 1024 300])
    hold all;
    for i=1:numel(x);
        bar_h=bar(x(i),n(i),0.3,'FaceColor',G2R(i,:));
    end
    xlabel('Co-reg error (mm)');
    ylabel('Count');
    title('Scalp - digitised head points co-reg error');
    set(gca,'FontSize',18,'FontWeight','bold')
    saveas(gcf,strcat(output_dir{1},'/RAC_head_points.png'));
    fprintf('Number of Points: %8.0f\nMean Displacement: %8.2fmm\nSTD of displacement: %6.2fmm\nMedian Displacement: %6.2fmm\nUpper Quartile: %11.2fmm\nMax Error: %16.2fmm\n',...
        numel(err),...
        mean(err),...
        std(err),...
        median(err),...
        quantile(err,[0.75]),...
        max(err))
end

fprintf('\n==================================================================\n')
fprintf('CHECK 9: Final check - was the scalp extraction any good?\n\n')
scalp_extraction_score=0;
if strcmp(coreg_method,'rhino')==1;
    mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);
    mesh_rhino_old=mesh_rhino;
    mesh_rhino = struct('faces',mesh_rhino.face,'vertices',mesh_rhino.vert);
    [a,b,c]=findTriMeshHoles(mesh_rhino.faces,mesh_rhino.vertices);
    if length(a)<=2;
        scalp_extraction_score=1;
        fprintf('Detected %1.0f holes. This is plausibly okay',round(length(a)));
    else
        fprintf('Detected %1.0f holes. We shouldn''t really have this many.',round(length(a)));
    end
else %spm
    mesh_spm=struct('faces',mesh.tess_scalp.face,'vertices',mesh.tess_scalp.vert);
    [a,b,c]=findTriMeshHoles(mesh_spm.faces,mesh_spm.vertices);
    if length(a)<=2;
        scalp_extraction_score=1;
        fprintf('Detected %1.0f holes. This is plausibly okay.\n',round(length(a)));
    else
        fprintf('Detected %1.0f holes. We shouldn''t really have this many.\n',round(length(a)));
    end
end

fprintf('\n==================================================================\n')
fprintf('All checks complete. Generating RAC report...\n')

%%% Make RAC report figure
make_RAC_report(D,fid_in_score,LR_vert_score,fid_delta_score,...
    brain_in_skull_score,skull_in_dewar_score,pol_point_score,...
    pol_points,custom_anat_score,err,scalp_extraction_score,output_dir);


%%% Make useful visualisation figure
if strcmp(coreg_method,'rhino')==1;
    visualise_co_reg(D,coreg_method,mesh_rhino,pol_points,scalp_extraction_score,a,output_dir)
else % we used SPM
    visualise_co_reg(D,coreg_method,mesh_spm,pol_points,scalp_extraction_score,a,output_dir)
end


%%% Write long report to file
total_score=0;
[~, fullname] = system('id -F');
fileID=fopen(strcat(output_dir{1},'/RAC_long_report.txt'),'w');
fprintf(fileID,'==================================================================');
fprintf(fileID,'\nRAC Report %s\n',string(datetime(now,'ConvertFrom','datenum')));
fprintf(fileID,'==================================================================\n');
fprintf(fileID,'Script run by user %s\n',fullname);
fprintf(fileID,'Original author: Ryan Timms, OHBA, University of Oxford\n');
fprintf(fileID,'ryan.timms@dtc.ox.ac.uk, @blobsonthebrain\n');
fprintf(fileID,'\nDataset information:\n');
fprintf(fileID,'Working on file %s\n',D.fullfile);
fprintf(fileID,'Co-reg method: %s\n',upper(coreg_method));
fprintf(fileID,'Scanner manufacturer: %s\n\n',upper(scanner_type));
fprintf(fileID,'-------------------\nChecklist Results:\n-------------------\n');
fprintf(fileID,'CHECK 1: Did the fiducials end up inside the dewar?\n');
if fid_in_score==3;
    fprintf(fileID,'All the fiducials ended up inside the MEG scanner\n');
    total_score=total_score+1;
elseif fid_in_score==0
    fprintf(fileID,'None of the fiducials ended up inside the MEG scanner!\n');
else
    fprintf(fileID,'The following fiducials ended up outside the dewar after co-reg: %s. Proceed with caution.\n',string(str_log));
end
fprintf(fileID,'\nCHECK 2: Are the L/RPA points aligned after co-reg?\n');
if LR_vert_score==1;
    total_score=total_score+1;
    fprintf(fileID,'The L/RPA points were roughly vertically aligned after co-reg.');
else
    fprintf(fileID,'The L/RPA points were *not* vertically aligned after co-reg.');
end
fprintf(fileID,' Vertical displacement: %fmm\n',vert_misalign);
fprintf(fileID,'\nCHECK 3: Did the fiducials end up within ±5mm of one another?\n');
if fid_delta_score==3;
    total_score=total_score+1;
    fprintf(fileID,'Yes\n%s: %fmm\n%s: %fmm\n%s: %fmm\n',string(fid_label_store(1)),delta_store(1),string(fid_label_store(2)),delta_store(2),string(fid_label_store(3)),delta_store(3));
else
    fprintf(fileID,'No\n%s: %fmm\n%s: %fmm\n%s: %fmm\n',string(fid_label_store(1)),delta_store(1),string(fid_label_store(2)),delta_store(2),string(fid_label_store(3)),delta_store(3));
end
fprintf(fileID,'--------------------------------\nConstituent fiducial errors:\n--------------------------------\n');
fprintf(fileID,'%s:\n%fmm\n%fmm\n%fmm\n',string(fid_label_store(1)),diff(1,1),diff(1,2),diff(1,3));
fprintf(fileID,'%s:\n%fmm\n%fmm\n%fmm\n',string(fid_label_store(2)),diff(2,1),diff(2,2),diff(2,3));
fprintf(fileID,'%s:\n%fmm\n%fmm\n%fmm\n',string(fid_label_store(3)),diff(3,1),diff(3,2),diff(3,3));
fprintf(fileID,'\nCHECK 4: Is the brain inside the skull?\n');
if brain_in_skull_score==2;
    total_score=total_score+1;
    fprintf(fileID,'The brain was found inside the skull\n');
elseif brain_in_skull_score==1;
    total_score=total_score+1;
    fprintf(fileID,'The vast majority of the brain was found inside the skull\n');
else
    fprintf(fileID,'The brain was found *outside* of the skull\n.');
end
fprintf(fileID,'\nCHECK 5: Was the skull in the dewar?\n');
if skull_in_dewar_score==2;
    total_score=total_score+1;
    fprintf(fileID,'Yes.\n');
elseif skull_in_dewar_score==1;
    fprintf(fileID,'Partly. This could be okay. Recommended you visualise the co-reg.\n');
else
    fprintf(fileID,'No. Make sure you visualise the co-reg!\n');
end
fprintf(fileID,'\nCHECK 6: Were any head points used in the co-reg?\n');
if isempty(pol_points)==1;
    fprintf(fileID,'No. Proceed with caution (unless you have a good reason for not having head points, e.g. you''re analysing headcast data)\n');
    fprintf(fileID,'Could not run CHECK 8 without these points\n');
else
    total_score=total_score+1;
    fprintf(fileID,'Yes! Summary statistics of digitised head points (CHECK 8) shown below:\n');
    fprintf(fileID,'Number of Points: %8.0f\nMean Displacement: %8.2fmm\nSTD of displacement: %6.2fmm\nMedian Displacement: %6.2fmm\nUpper Quartile: %11.2fmm\nMax Error: %16.2fmm\n',...
        numel(err),...
        mean(err),...
        std(err),...
        median(err),...
        quantile(err,[0.75]),...
        max(err));
end
fprintf(fileID,'\nCHECK 7: Was a custom MRI used to do the co-reg?\n');
if custom_anat_score==1;
    total_score=total_score+1;
    fprintf(fileID,'Yes. Good.\n');
else
    fprintf(fileID,'No. Why not? Be *very* careful in interpreting results if you plan to do source localisation.\n');
end
fprintf(fileID,'\nCHECK 9: Final check - was the scalp extraction any good?\n');
if scalp_extraction_score==0;
    fprintf(fileID,'This looks somewhat sub-optimal. We detected %2.0f holes in the scalp extraction patch. Take the co-reg with a pinch of salt.\n',round(length(a)));
else
    total_score=total_score+1;
    fprintf(fileID,'Didn''t detect any holes in the output of the scalp extraction.\n');
end
fprintf(fileID,'\nConclusion:\n');
if total_score==8;
    fprintf(fileID,'From the metrics computed here, it looks like the result of the co-reg was good.\n However, we always strongly advocate that you visualise your co-reg at least once.\n Please do so now if you haven''t already!');
elseif total_score==0;
    fprintf(fileID,'This co-reg was probably a disaster. You *must* investigate what has happened!');
elseif total_score<8
    fprintf(fileID,'Hard to conclude whether or not this co-reg was a complete success without you looking at it. This co-reg passed %2.0f of the 8/9 checks.\n',round(total_score));
end
fclose(fileID);

%%% Write short report to file
fileID=fopen(strcat(output_dir{1},'/RAC_short_report.txt'),'w');
fprintf(fileID,'RAC report for %s\n',D.fullfile);
fprintf(fileID,'%s\n',string(datetime(now,'ConvertFrom','datenum')));
fprintf(fileID,'Co-reg method: %s\n',upper(coreg_method));
fprintf(fileID,'Scanner manufacturer: %s\n\n',upper(scanner_type));
fprintf(fileID,'-------------------\n-------------------\nChecklist Results:\n-------------------\n-------------------\n');
fprintf(fileID,'fid_in_score: %s\n',string(fid_in_score));
fprintf(fileID,'LR_vert_score: %s\n',string(LR_vert_score));
fprintf(fileID,'Vertical displacement: %fmm\n',vert_misalign);
fprintf(fileID,'fid_delta_score: %s\n',string(fid_delta_score));
fprintf(fileID,'--------------------------------\nRMS fiducial errors:\n--------------------------------\n');
fprintf(fileID,'%s: %fmm\n%s: %fmm\n%s: %fmm\n',string(fid_label_store(1)),delta_store(1),string(fid_label_store(2)),delta_store(2),string(fid_label_store(3)),delta_store(3));
fprintf(fileID,'--------------------------------\nConstituent fiducial errors:\n--------------------------------\n');
fprintf(fileID,'%s:\n%fmm\n%fmm\n%fmm\n',string(fid_label_store(1)),diff(1,1),diff(1,2),diff(1,3));
fprintf(fileID,'%s:\n%fmm\n%fmm\n%fmm\n',string(fid_label_store(2)),diff(2,1),diff(2,2),diff(2,3));
fprintf(fileID,'%s:\n%fmm\n%fmm\n%fmm\n',string(fid_label_store(3)),diff(3,1),diff(3,2),diff(3,3));
fprintf(fileID,'brain_in_skull_score: %s\n',string(brain_in_skull_score));
fprintf(fileID,'skull_in_dewar_score: %s\n',string(skull_in_dewar_score));
fprintf(fileID,'pol_points_used: %s\n',string(abs(1-isempty(pol_points))));
fprintf(fileID,'custom_anat_score: %s\n',string(custom_anat_score));
fprintf(fileID,'scalp_extraction_score: %s\n',string(scalp_extraction_score));
fprintf(fileID,'Number of holes: %s',string(round(length(a))));
if sum(isnan(err))==0;
    fprintf(fileID,'\nSummary statistics of digitised head points (CHECK 8) shown below:\n');
    fprintf(fileID,'Number of Points: %8.0f\nMean Displacement: %8.2fmm\nSTD of displacement: %6.2fmm\nMedian Displacement: %6.2fmm\nUpper Quartile: %11.2fmm\nMax Error: %16.2fmm\n',...
        numel(err),...
        mean(err),...
        std(err),...
        median(err),...
        quantile(err,[0.75]),...
        max(err));
else
    fprintf('\n');
end


fprintf(fileID,'\n-------------------\nKey/formatting:\n-------------------\nfid_in_score - number of fiducials inside the dewar after co-reg. Ranges from 0 to 3.\nLR_vert_score - 0 or 1, corresponding to whether fiducials are (1) or are not (0) aligned.\nVertical displacement - vertical misalignment of fiducials after co-reg.\n');
fprintf(fileID,'fid_delta_score - Did the fiducials end up within ±5mm of one another?. Ranges from 0 (all have a delta >5mm) to 3 (all displacement vectors <5mm).');
fprintf(fileID,'\nLPA/RPA/Nasion - Euclidian distance between fiducials after the co-reg. Note that these can appear in different order and have slightly different naming conventions (e.g. "Nasion" vs "NAS").');
fprintf(fileID,'\nbrain_in_skull_score - Was the brain found inside the skull? [0/1/2] for no/most of it/yes]');
fprintf(fileID,'\nskull_in_dewar_score - 0, 1 or 2, corresponding to skull being fully outside, partly outside or completely inside the dewar, respectively.');
fprintf(fileID,'\npol_points_used - 0 (no points used) or 1 (head points used). If equal to 1, statistics are printed at the end of the log file.');
fprintf(fileID,'\ncustom_anat_score - subject anatomy not used/used [0/1].');
fprintf(fileID,'\nscalp_extraction_score - 0 or 1, corresponding to holes being found in the scalp extraction.');
fclose(fileID);
fprintf('\nDone. Reports and Figures saved to disk.\nCiao.\n');
close all; % goodbye
end

function in = inhull(testpts,xyz,tess,tol)
% inhull: tests if a set of points are inside a convex hull
% usage: in = inhull(testpts,xyz)
% usage: in = inhull(testpts,xyz,tess)
% usage: in = inhull(testpts,xyz,tess,tol)
%
% arguments: (input)
%  testpts - nxp array to test, n data points, in p dimensions
%       If you have many points to test, it is most efficient to
%       call this function once with the entire set.
%
%  xyz - mxp array of vertices of the convex hull, as used by
%       convhulln.
%
%  tess - tessellation (or triangulation) generated by convhulln
%       If tess is left empty or not supplied, then it will be
%       generated.
%
%  tol - (OPTIONAL) tolerance on the tests for inclusion in the
%       convex hull. You can think of tol as the distance a point
%       may possibly lie outside the hull, and still be perceived
%       as on the surface of the hull. Because of numerical slop
%       nothing can ever be done exactly here. I might guess a
%       semi-intelligent value of tol to be
%
%         tol = 1.e-13*mean(abs(xyz(:)))
%
%       In higher dimensions, the numerical issues of floating
%       point arithmetic will probably suggest a larger value
%       of tol.
%
%       DEFAULT: tol = 0
%
% arguments: (output)
%  in  - nx1 logical vector
%        in(i) == 1 --> the i'th point was inside the convex hull.
%
% Example usage: The first point should be inside, the second out
%
%  xy = randn(20,2);
%  tess = convhulln(xy);
%  testpoints = [ 0 0; 10 10];
%  in = inhull(testpoints,xy,tess)
%
% in =
%      1
%      0
%
% A non-zero count of the number of degenerate simplexes in the hull
% will generate a warning (in 4 or more dimensions.) This warning
% may be disabled off with the command:
%
%   warning('off','inhull:degeneracy')
%
% See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 10/26/06
% get array sizes
% m points, p dimensions
p = size(xyz,2);
[n,c] = size(testpts);
if p ~= c
    error 'testpts and xyz must have the same number of columns'
end
if p < 2
    error 'Points must lie in at least a 2-d space.'
end
% was the convex hull supplied?
if (nargin<3) || isempty(tess)
    tess = convhulln(xyz);
end
[nt,c] = size(tess);
if c ~= p
    error 'tess array is incompatible with a dimension p space'
end
% was tol supplied?
if (nargin<4) || isempty(tol)
    tol = 0;
end
% build normal vectors
switch p
    case 2
        % really simple for 2-d
        nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * [0 1;-1 0];
        
        % Any degenerate edges?
        del = sqrt(sum(nrmls.^2,2));
        degenflag = (del<(max(del)*10*eps));
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate edges identified in the convex hull'])
            
            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
    case 3
        % use vectorized cross product for 3-d
        ab = xyz(tess(:,1),:) - xyz(tess(:,2),:);
        ac = xyz(tess(:,1),:) - xyz(tess(:,3),:);
        nrmls = cross(ab,ac,2);
        degenflag = false(nt,1);
    otherwise
        % slightly more work in higher dimensions,
        nrmls = zeros(nt,p);
        degenflag = false(nt,1);
        for i = 1:nt
            % just in case of a degeneracy
            % Note that bsxfun COULD be used in this line, but I have chosen to
            % not do so to maintain compatibility. This code is still used by
            % users of older releases.
            %  nullsp = null(bsxfun(@minus,xyz(tess(i,2:end),:),xyz(tess(i,1),:)))';
            nullsp = null(xyz(tess(i,2:end),:) - repmat(xyz(tess(i,1),:),p-1,1))';
            if size(nullsp,1)>1
                degenflag(i) = true;
                nrmls(i,:) = NaN;
            else
                nrmls(i,:) = nullsp;
            end
        end
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate simplexes identified in the convex hull'])
            
            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
end
% scale normal vectors to unit length
nrmllen = sqrt(sum(nrmls.^2,2));
% again, bsxfun COULD be employed here...
%  nrmls = bsxfun(@times,nrmls,1./nrmllen);
nrmls = nrmls.*repmat(1./nrmllen,1,p);
% center point in the hull
center = mean(xyz,1);
% any point in the plane of each simplex in the convex hull
a = xyz(tess(~degenflag,1),:);
% ensure the normals are pointing inwards
% this line too could employ bsxfun...
%  dp = sum(bsxfun(@minus,center,a).*nrmls,2);
dp = sum((repmat(center,nt,1) - a).*nrmls,2);
k = dp<0;
nrmls(k,:) = -nrmls(k,:);
% We want to test if:  dot((x - a),N) >= 0
% If so for all faces of the hull, then x is inside
% the hull. Change this to dot(x,N) >= dot(a,N)
aN = sum(nrmls.*a,2);
% test, be careful in case there are many points
in = false(n,1);
% if n is too large, we need to worry about the
% dot product grabbing huge chunks of memory.
memblock = 1e6;
blocks = max(1,floor(n/(memblock/nt)));
aNr = repmat(aN,1,length(1:blocks:n));
for i = 1:blocks
    j = i:blocks:n;
    if size(aNr,2) ~= length(j),
        aNr = repmat(aN,1,length(j));
    end
    in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
end
end


function [holeCellArray,bounding_triangles,holeLengths] = findTriMeshHoles(faces,vertices)
% Finds holes in a triangular mesh
% Note: Does not work if a hole shares more than one vertex with other holes
% Input:
% faces = M x 3
% vertices = N x 3 (optional if you want the hole lengths)
% Output:
% holeCellArray = P x 1 cell array containing a list of holes, which are
% traced in consecutive order (list of scalar indices)
% bounding_triangles = Q x 3 list of faces that contain a bounding edge (does
% not contain triangles that only has a single bounding vertex)
% holeLengths = P x 1 vector containing the perimeter of each hole
if nargout > 2 && nargin < 2
    error('Please input vertices if you want the hole lengths')
end
% sort each row
faces = sort(faces,2);
% face error checking (ensure no faces contain an index more than once)
sel = faces(:,1) == faces(:,2) | faces(:,2) == faces(:,3) | faces(:,1) == faces(:,3);
if sum(sel) > 0
    warning('Triangle refers to a vertex more than once... Consider fixing your data...')
    faces = faces(~sel,:);
end
% get edges from faces
edge1 = faces(:,1:2);
edge2 = faces(:,[1,3]);
edge3 = faces(:,2:3);
edges = vertcat(edge1,edge2,edge3);
% find edges that have no replicates (find boundary edges)
edges = sortrows(edges);
vec = diff(edges); % diff btw rows
vec = any(vec,2); % unique rows
vec1 = [1;vec]; % doesn't match with previous row
vec2 = [vec;1]; % doesn't match with next row
vec = vec1&vec2; % doesn't match previous and next row
bounding_edges = edges(vec,:); % unique edges
% bounding triangles - triangles that contain a bounding edge
% Note: you can modify to get triangles that contain a bounding vertex
match1 = ismember(edge1,bounding_edges,'rows');
match2 = ismember(edge2,bounding_edges,'rows');
match3 = ismember(edge3,bounding_edges,'rows');
bounding_triangles = faces(match1|match2|match3,:);
if nargin > 1 && nargout > 1
    % length of each edge
    edgeVertex1 = vertices(bounding_edges(:,1),:);
    edgeVertex2 = vertices(bounding_edges(:,2),:);
    edge_lengths = sqrt(sum((edgeVertex1-edgeVertex2).^2,2));
    [holeCellArray,holeLengths]=getTraces(bounding_edges,edge_lengths);
else
    holeCellArray=getTraces(bounding_edges);
end
end
function [tracesDB,lengthDB,completeDB]=getTraces(bounding_edges,edge_lengths)
% Connects the unordered edges to form holes that list the edges in
% consecutive order.
% The number of formed holes from the edges determines the size of tracesDB
% INPUT
% bounding_edges = N x 2 list of connected border vertices
% edge_lengths   = N x 1 the corresponding length of the connected vertices
% OUTPUT
%   tracesDB = cell array containing a list of holes, which are traced in
%               consecutive order
%   lengthDB = 1-D vector containing the perimeter of each hole that
%                   corresponds to tracesDB
%   completeDB = 1-D vector that indicates if corresponding hole in tracesDB is completely traced
%                 (technically should be all 1's) - for testing purposes
if nargin < 2
    edge_lengths = zeros(size(bounding_edges,1),1);
end
tracesDB = {};
lengthDB = [];
completeDB = [];
trace12 = [];
idx = 1;
while size(bounding_edges,1) > 0
    trace = bounding_edges(1,:)'; % add first trace
    bounding_edges = bounding_edges(2:end,:);
    total_length = edge_lengths(1); % add first trace length
    edge_lengths = edge_lengths(2:end);
    flip_flag = 0;
    complete_flag = 1;
    while trace(1) ~= trace(end)
        [row_idx,col_idx] = find(bounding_edges == trace(end));
        col_idx = -col_idx+3;
        % if more than one bounding edge connects to the end of the
        % trace then trace in the other direction (occurs when only one
        % triangle separates two holes so both holes share a single vertex)
        if length(row_idx) > 1 || length(row_idx) < 1
            if flip_flag == 0
                trace = flip(trace);
                flip_flag = 1;
                continue;
            else
                complete_flag = 0;
                %             disp('Incomplete hole tracing... May need to make more conditions')
                break;
            end
        end
        trace = [trace;bounding_edges(row_idx,col_idx)];
        total_length = total_length+edge_lengths(row_idx);
        bounding_edges(row_idx,:) = [];
        edge_lengths(row_idx) = [];
    end
    tracesDB = [tracesDB;{trace}];
    lengthDB = [lengthDB; total_length];
    completeDB = [completeDB;complete_flag];
    trace12 = [trace12;trace(1),trace(end)];
    idx = idx+1;
end
if sum(~completeDB) > 0
    warning('Not all holes are completely traced due to holes sharing more than one vertex with other holes')
end
end

function RT_rhino_display(D,hf)
% Display of the coregistered RHINO meshes and polhemus/sensor locations.
%
%
% RHINO_DISPLAY(D)    - displays the coregistration for MEEG object D in a
%                       new figure
%
% RHINO_DISPLAY(D,hf) - displays the coregistration for MEEG object D in an
%                       existing figure with handle hf
%
% Adam Baker 2015
% Small modifications by Ryan Timms, 2020, to make compatible with OSL_RAC



hold all
D = montage(D,'switch');

mesh       = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
mesh_rhino = spm_eeg_inv_transform_mesh(D.inv{1}.datareg(1).fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh.tess_rhino);

% Reduce mesh size for faster plotting
mesh_rhino = struct('faces',mesh_rhino.face,'vertices',mesh_rhino.vert);
%reduce_factor = 1e3/size(mesh_rhino.vertices,1);
%mesh_rhino = reducepatch(mesh_rhino,reduce_factor);

headshape_polhemus = D.inv{1}.datareg(1).fid_eeg.pnt;
fid_polhemus = D.inv{1}.datareg(1).fid_eeg.fid.pnt;
fid_mri = D.inv{1}.datareg(1).fid_mri.fid.pnt;

patch(mesh_rhino,'FaceColor',[238,206,179]./255,'EdgeColor','none','FaceAlpha',0.6);
patch(struct('faces',mesh.tess_ctx.face,'vertices',mesh.tess_ctx.vert),'FaceColor',[1 0.7 0.7],'EdgeColor','none');
%patch(struct('faces',mesh.tess_scalp.face,'vertices',mesh.tess_scalp.vert),'EdgeColor','g','FaceColor','none','Parent',ha);



axis('image','off')
material shiny
lighting gouraud
hl = camlight('headlight');
% set(hf,'WindowButtonMotionFcn',@(~,~) camlight(hl,'headlight'),'Interruptible','off','busyaction','cancel');


% Headshape locations
%--------------------------------------------------------------------------
if ~isempty(headshape_polhemus)
    %h_hsp   = plot3(headshape_polhemus(:,1),headshape_polhemus(:,2),headshape_polhemus(:,3),'dm');
    %set(h_hsp,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','r');
    
    err = rhino_cperror(mesh_rhino.vertices,headshape_polhemus);
    err(err>5) = 5;
    h_hsp = scatter3(headshape_polhemus(:,1),headshape_polhemus(:,2),headshape_polhemus(:,3),50,err,'filled');
    colormap(jet)
end

% Sensors (coreg.)
%--------------------------------------------------------------------------
h_sens  = plot3(D.sensors('MEG').chanpos(:,1),D.sensors('MEG').chanpos(:,2),D.sensors('MEG').chanpos(:,3),'og');
set(h_sens,'MarkerFaceColor','g','MarkerSize', 12,'MarkerEdgeColor','k');

% EEG fiducials or MEG coils (coreg.)
%--------------------------------------------------------------------------
h_fid   = plot3(fid_polhemus(:,1),fid_polhemus(:,2),fid_polhemus(:,3),'oc');
set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');

% MRI fiducials
%--------------------------------------------------------------------------
h_fidmr = plot3(fid_mri(:,1),fid_mri(:,2),fid_mri(:,3),'dm');
set(h_fidmr,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');

material dull % Improve visibility of the brain by reducing reflection glare
axis off


end

function make_RAC_report(D,fid_in_score,LR_vert_score,fid_delta_score,...
    brain_in_skull_score,skull_in_dewar_score,pol_point_score,...
    pol_points,custom_anat_score,err,scalp_extraction_score,output_dir);
%%% Make RAC report.
figure('Position',[10 10 1600 1024])
subplot(10,4,[1:4]);
text(0,10,char(strcat('RAC Report for dataset',{' '},string(D.fname))),'FontSize',24,'Interpreter','none','FontWeight','bold')
ylim([0 20])
axis off

subplot(10,4,[1:3]+4);
% CHECK 1: Are the fiducials inside the scanner?
if fid_in_score==3;
    text(0,10,'All the fiducials ended up inside the MEG scanner','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[4]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
elseif fid_in_score==0
    text(0,10,'None of the fiducials ended up inside the MEG scanner!','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[4]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
else
    text(0,10,'Some of the fiducials ended up outside of the MEG scanner. Proceed with caution','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[4]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[255, 102, 0]./255,'EdgeColor','k','LineWidth',2)
    axis equal
    axis off
end

% Check 2: are the L/RPA points aligned after co-reg?
if LR_vert_score==1;
    subplot(10,4,[5:7]+4);
    text(0,10,'The L/RPA points were roughly vertically aligned after co-reg','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[8]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
else
    subplot(10,4,[5:7]+4);
    text(0,10,'The L/RPA points were *not* vertically aligned after co-reg','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[8]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
end

if fid_delta_score==3;
    subplot(10,4,[9:11]+4);
    text(0,10,'The MEG/MR fiducials ended up being within ±5mm of one another after co-reg','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[12]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
else
    subplot(10,4,[9:11]+4);
    text(0,10,'The MEG/MR fiducials were *not* within ±5mm of one another after co-reg','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[12]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
end

if brain_in_skull_score==2;
    subplot(10,4,[13:15]+4);
    text(0,10,'The brain was found to be inside the skull','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[16]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
elseif brain_in_skull_score==1;
    subplot(10,4,[13:15]+4);
    text(0,10,'The vast majority of the brain was found in the skull. Some was outside.','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[16]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[255, 102, 0]./255,'EdgeColor','k','LineWidth',2)
    axis equal
    axis off
     
else
    text(0,10,'The brain was *not* found in the skull','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[16]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
end

% Check 5: skull in dewar
if skull_in_dewar_score==2;
    subplot(10,4,[17:19]+4);
    text(0,10,'The skull was completely in the dewar','FontSize',18)
    ylim([0 20])
    axis off
    
    
    subplot(10,4,[20]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
elseif skull_in_dewar_score==0
    subplot(10,4,[17:19]+4);
    text(0,10,'The skull was not in the dewar at all','FontSize',18)
    ylim([0 20])
    axis off
    
    
    subplot(10,4,[20]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
else
    subplot(10,4,[17:19]+4);
    text(0,10,'At least 80 percent of skull points were in the dewar. This is probably okay.','FontSize',18)
    ylim([0 20])
    axis off
    
    
    subplot(10,4,[20]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[255, 102, 0]./255,'EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
end

if pol_point_score==0;
    subplot(10,4,[21:23]+4);
    text(0,10,'No head points used for co-reg!','FontSize',18)
    ylim([0 20])
    axis off
    
    
    subplot(10,4,[24]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
else
    subplot(10,4,[21:23]+4);
    text(0,10,'Head points were used for co-reg!','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[24]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
end

if custom_anat_score==0;
    subplot(10,4,[25:27]+4);
    text(0,10,'No personalised MRI used','FontSize',18)
    ylim([0 20])
    axis off
    
    
    subplot(10,4,[28]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','r','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
else
    subplot(10,4,[25:27]+4);
    text(0,10,'A personalised structural MR was used for co-reg','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[28]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
end

if isempty(pol_points)==1;
    subplot(10,4,[29:31]+4);
    text(0,10,'Polhemus statistics N/A','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[32]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','w','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
elseif median(err)<5 && max(err)<10;
    subplot(10,4,[29:31]+4);
    text(0,10,['Median scalp-head-point displacement is',num2str(median(err)),'mm'],'FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[32]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
    
elseif median(err)<5 && max(err)>10;
    subplot(10,4,[29:31]+4);
    text(0,10,['Median scalp-head-point displacement is ',num2str(median(err),3),'mm, but also detected large outlier(s) >10mm'],'FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[32]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[255, 102, 0]./255,'EdgeColor','k','LineWidth',2)
    axis equal
    axis off
else
    subplot(10,4,[29:31]+4);
    text(0,10,['Median scalp-head-point displacement is ',num2str(median(err),3),'mm'],'FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[32]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[255, 102, 0]./255,'EdgeColor','k','LineWidth',2)
    axis equal
    axis off
end

if scalp_extraction_score==0;
    subplot(10,4,[33:35]+4);
    text(0,10,'Detected more than 2 holes in scalp extracted surface. Proceed with caution','FontSize',18)
    ylim([0 20])
    axis off
    
    
    subplot(10,4,[36]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor',[255, 102, 0]./255,'EdgeColor','k','LineWidth',2)
    axis equal
    axis off
else
    subplot(10,4,[33:35]+4);
    text(0,10,'Scalp extraction was clean','FontSize',18)
    ylim([0 20])
    axis off
    
    subplot(10,4,[36]+4);
    rectangle('Position',[0 0 1 1],'Curvature',[1 1],'FaceColor','g','EdgeColor','k','LineWidth',2)
    axis equal
    axis off
    
end
saveas(gcf,strcat(output_dir{1},'/RAC_traffic_report.png'));
end

function visualise_co_reg(D,coreg_method,mesh,pol_points,scalp_extraction_score,a,output_dir)

% Inputs:
%   D object
%   Co-reg method
%   mesh: rhino or spm, depending on co-reg method
%   pol_points: digitised head points
%   scalp_extraction_score: flag as to whether or not BETSURF was good.
%   a: cell-array containing hole locations

if strcmp(coreg_method,'rhino')==1;
    mesh_rhino=mesh;
else
    mesh_spm=mesh;
end

figure('Position',[0 0 2048 1024])
subplot(4,3,[1:3])
text(0,1,char(strcat('Co-reg result for dataset',{' '},string(D.fname))),'FontSize',25,'Interpreter','none','FontWeight','bold')
hold all
xlim([0 1]);
ylim([0 2]);
axis off

% Right view
subplot(4,3,[4,7,10]);
if strcmp(coreg_method,'rhino')==1;
    RT_rhino_display(D)
    view([90 0]);
else
    spm_eeg_inv_checkdatareg_3Donly(D);
    view([0 0]);
end

% Also display the holes on this right-hand view figure
if scalp_extraction_score==0;
    if strcmp(coreg_method,'rhino')==1;
        for i=1:length(a)
            plot3(mesh_rhino.vertices(a{i},1),mesh_rhino.vertices(a{i},2),mesh_rhino.vertices(a{i},3),'r-','LineWidth',2)
        end
    else
        plot3(mesh_spm.vertices(a{i},1),mesh_spm.vertices(a{i},2),mesh_spm.vertices(a{i},3),'r-','LineWidth',2)
    end
end

zoom(0.5)
title('Right','FontSize',18);

% Top view
subplot(4,3,[5,8,11]);
if strcmp(coreg_method,'rhino')==1;
    RT_rhino_display(D)
else
    spm_eeg_inv_checkdatareg_3Donly(D);
    view([-90 90]);
    
end
zoom(0.5);
title('Top','FontSize',18);

% Front view
subplot(4,3,[6,9,12]);
if strcmp(coreg_method,'rhino')==1;
    RT_rhino_display(D)
    view([180 0]);
else
    spm_eeg_inv_checkdatareg_3Donly(D);
    view([90 0]);
    zoom(0.5);
end
zoom(0.9)
title('Front','FontSize',18);

% Legends are bloody awkward in MATLAB:
h=get(gca,'Children'); % grab all the axes handles at once
if strcmp(coreg_method,'rhino')==1;
    legendstr={'MR Fiducials','MEG Fiducials','MEG Sensors','Digitised Head Points','ee'};
    if isempty(pol_points)~=0;
        lgd=legend(h([1 2 3 4]),legendstr{[1 2 3 4]});
    else
        lgd=legend(h([1 2 3]),legendstr{[1 2 3]});
    end
else
    legendstr={'MR Fiducials','MEG Fiducials','MEG Sensors','Digitised Head Points','ee'};
    if isempty(pol_points)~=1;
        lgd=legend(h([1 2 3 4]),legendstr{[1 2 3 4]});
    else
        lgd=legend(h([1 2 3]),legendstr{[1 2 3]});
    end
end
lgd.Position=[0.8072,0.822,0.1189,0.1074];
lgd.FontSize=20;
title(lgd,'Key','FontSize',20)
hold off
if isempty(pol_points)~=1;
    c=colorbar('Location','North');
    c.FontSize=20;
    c.TickDirection='in';
    c.Position=[0.364510869565217,0.012845331720105*5,0.310489130434783,0.035438596491228];
    ylabel(c,'Distance (mm)');
end
try
    saveas(gcf,strcat(output_dir{1},'/RAC_coreg_vis.png'));
catch
    warning('Unable to save RAC report Figure to disk.');
end
try
    saveas(gcf,strcat(output_dir{1},'/RAC_coreg_vis.fig'));
catch
    warning('Unable to save RAC report Figure to disk.');
end
end

function [max_x,min_x,max_y,min_y,max_z,min_z]=dewar_extremities(gradiometer_locations);
% Approximate the edges of the dewar by finding the most extreme sensor
% locations. N.b. that this assumes that we have removed reference coils
% from 'gradiometer_locations' in the case of CTF data

[where,~]=max(gradiometer_locations(:,1)); %x
max_x=where;

[where,~]=min(gradiometer_locations(:,1));
min_x=where;

[where,~]=max(gradiometer_locations(:,2)); %y
max_y=where;

[where,~]=min(gradiometer_locations(:,2));
min_y=where;

[where,~]=max(gradiometer_locations(:,3)); %z
max_z=where;

[where,~]=min(gradiometer_locations(:,3));
min_z=where;
end