function [D,parcellation,assignments,mni_coords] = osl_apply_parcellation(S)
% Computes node time series from a source space MEEG object using a
% parcellation
%
% [D,parcellation,assignments,mni_coords] = osl_apply_parcellation(S)
%
% INPUTS:
%
% S.D                   - source space MEEG object
% S.parcellation        - .nii or .nii.gz file at the same gridstep as S.D
%                         OR a matrix of dimensions voxels x parcels
%                         with the same number of voxels as S.D
% S.orthogonalisation   - Orthogonalisation protocol to apply:
%                           ['none','symmetric','closest','householder']
% S.method              - method for reconstructing parcel time course:
%                           ['PCA','mean','peakVoxel']
% S.prefix              - filename prefix for new MEEG object (default 'p')
% S.hcp_sourcemodel3d   - if it is HCP data, then need to pass in the
%                           sourcemodel3d
% S.trialwise           - set to 1 to compute orth separately for each
%                         trial, default is 0
%
% OUTPUTS:
%
% D             - new MEEG object containing the parcel time courses
% parcellation  - the parcellation weights for each voxel
% assignments   - parcel assignments for each voxel from a winner takes all
%                 voting
% mni_coords    - mni_coords of centres of parcels

% Check SPM File Specification:
try
    if isa(S.D,'meeg')
        D = S.D;
    else
        S.D = char(S.D);
        [pathstr,filestr] = fileparts(S.D);
        S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
        D = spm_eeg_load(S.D);
    end;
    D.check;
catch ME
    disp(ME.message);
    error(['SPM file specification not recognised or incorrect']);
end

try
    D = montage(D,'switch',2);
catch ME
    display(ME.message);
    display('Continuing with existing montage');
end
    
S.prefix            = ft_getopt(S,'prefix','p');
S.orthogonalisation = ft_getopt(S,'orthogonalisation','none');
S.method            = ft_getopt(S,'method','PCA');

try 
    maskfname=S.maskfname;
catch ME
    display(ME.message);
    maskfname=[osldir '/std_masks/MNI152_T1_' num2str(getmasksize(D.nchannels)) 'mm_brain.nii.gz'];                
end

if isfield(S,'hcp_sourcemodel3d') && ~isempty(S.hcp_sourcemodel3d)
    % HCP data
    parcellation = hcp_mni2hcp_giles(S.hcp_sourcemodel3d,S.parcellation,S.hcp_mask_fname_out);
else
    switch class(S.parcellation)
        case {'char','cell'}
            try
                stdbrain = nii.load(maskfname);
                parc=nii.load(S.parcellation);
                parcellation = vols2matrix(parc,stdbrain); %nVoxels x nSamples
            catch
                error('Make sure the parcellation file and the data are valid and compatible, including having the same spatial resolution.');
            end

            % compute mni_coords as centres of gravity of parcels
            mni_coords=osl_mniparcellation2mnicoords(S.parcellation);

        case {'single','double','logical'}
            parcellation = S.parcellation;
        otherwise
            error('Unrecognized parcellation');
    end
end

try 
    innovations_mar_order=S.innovations_mar_order;
catch
    innovations_mar_order=12;
end

if size(parcellation,2) == 1
    % Currently nvoxels x 1 with a index at each voxel indicating
    % parcel membership. Instead, needs to be nvoxels x nparcels binary labelling
    num_parcels = max(parcellation);
    parcellation_new = zeros(size(parcellation,1),num_parcels);
    for voxel = 1:size(parcellation,1)
        if parcellation(voxel) > 0
            parcellation_new(voxel,parcellation(voxel)) = 1;
        end
    end
    parcellation = parcellation_new;
    clear parcelflag_new;
end

% create parcellation based on winner takes all voting
% (this will be used for outputting results in nii format)
assignments = false(size(parcellation));
for voxel = 1:size(parcellation,1)
    [a,b] = max(parcellation(voxel,:));
    if a > 0
        assignments(voxel,b) = 1;
    end
end


if D.ntrials == 1 % can just pass in the MEEG object
    voxeldata = D(:,:,:);
    goodsamples = good_samples(D); % all checks if any channels are bad then ignore those samples
    % Get time-coursess
    nodedata = ROInets.get_node_tcs(voxeldata(:,goodsamples), parcellation, S.method);
    if ~strcmp(S.orthogonalisation,'none')
        
        if ~strcmp(S.orthogonalisation, 'innovations_mar')
            nodedata = ROInets.remove_source_leakage(nodedata,S.orthogonalisation);
        else            
            nodedata = leakcorr(nodedata',size(nodedata,2),innovations_mar_order)';
        end
    end

    data = zeros(size(nodedata,1),length(goodsamples));
    data(:,goodsamples) = nodedata;

    data = reshape(data,[size(data,1),length(goodsamples),D.ntrials]);
                
elseif isa(D,'meeg') % work with D object in get_node_tcs
    %goodsamples = find(good_samples(D));
    nodedata = ROInets.get_node_tcs(D,parcellation,S.method);
    nodedata = reshape(nodedata(:,:,:),size(nodedata,1),[]);
    if ~strcmp(S.orthogonalisation,'none')
        if ~strcmp(S.orthogonalisation, 'innovations_mar')
            nodedata = ROInets.remove_source_leakage(nodedata,S.orthogonalisation);
        else            
            nodedata = leakcorr(nodedata',size(nodedata,2),innovations_mar_order)';
        end
                  
        
    end
    data = nodedata;
    data = reshape(data,size(data,1),size(D,2),size(D,3));

else % reshape the data first (or fix get_node_tcs to work with trialwise MEEG data)
    goodsamples = good_samples(D);
    nodedata = zeros(size(parcellation,2),D.nsamples,D.ntrials);
    msg = '';

    for idx = 1:D.ntrials
        % We're ignoring trials with bad samples as they might be rank
        % deficient, just leave them empty.
        fprintf(repmat(char(8),1,length(msg)));
        msg = sprintf('Orthoganalising trial: %d of %d',idx,D.ntrials);
        fprintf(msg);
        if goodsamples(1,1,idx) == 1
            try
                voxeldata = D(:,:,idx);
                nodedata(:,:,idx) = ROInets.get_node_tcs(voxeldata, parcellation, S.method,0);
                if ~strcmp(S.orthogonalisation,'none')
                    if ~strcmp(S.orthogonalisation, 'innovations_mar')
                        nodedata = ROInets.remove_source_leakage(nodedata,S.orthogonalisation);
                    else            
                        nodedata(:,:,idx) = leakcorr(nodedata(:,:,idx)',size(nodedata,2),innovations_mar_order)';
                    end
                end
            catch,
                % This trial is probably low rank, ignore it
                goodsamples(:,:,idx) = 0;
                disp('Skipping low rank trial');
            end

        end

        data = nodedata;
        fprintf('\n');

    end
end

clear voxeldata_concat nodedata_concat;

% Save data to new MEEG object
outfile = prefix(fullfile(D.path,D.fname),S.prefix);
Dnode = clone(montage(D,'switch',0),outfile,[size(data,1),D.nsamples,D.ntrials]);
parcel_chantype='VE';
Dnode = chantype(Dnode,1:Dnode.nchannels,parcel_chantype);
Dnode(:,:,:) = data;

% copy badtrials
badtrials=D.badtrials;
if ~isempty(badtrials)
    Dnode = Dnode.badtrials(1:length(badtrials),badtrials);
end

% copy events, but need to change ev.value to be 'VE' for all artefact
% events, so that bad segments get passed through
ev = D.events;
for ee=1:length(ev)
    if strncmp(ev(ee).type,'artefact',8)
        ev(ee).value=parcel_chantype;
    end
end
Dnode = Dnode.events(1,ev);

Dnode.save;

D = Dnode; % For output

D.parcellation.weights=parcellation;
D.parcellation.assignments=assignments;
D.parcellation.S=S;
D.parcellation.mni_coords=mni_coords;

D.save;
end
