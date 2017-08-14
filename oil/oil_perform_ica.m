%% osl_perform_ica.m
%
%  res=osl_perform_ica(S)
%  
%  Input structure S contains:
%     REQUIRED: ica_concat_path: the full path to the concatenated data
%     (usually contructed using osl_concat_subs).
%               oil - the oil structure that defines the analysis
%               parameters.
%     OPTIONAL: num_ics        :  the number of independent components sought(default 25).   
%               last_eig       :  the number ofprinciple components used (default equal to num_ics).
%               icasso_its     :  number of ICASSO iterations (default is 0).
%               nonlinearity   :  fastICA nonlinearity function (default is 'tanh').               
%               temp_or_spat   :  either 'temporal' or 'spatial' ICA 
%               use_gm_mask    :  limits ICA to grey matter voxels only (default is 0).
%               gridstep       :  spatial resolution of ICA inout data.
%               normalise_vox  :  normalise individual subject voxels
%                                 (default = 0). (See
%                                 http://dx.doi.org/10.1016/j.neuroimage.2012.11.011 
%                                 for a discussion of why not to 
%                                 variance-normalise). 
%               subj_ind       :  individual subject separation indices (needed for voxel normalisation).
%
%  Output structure res.results contains:      
%               mixing_matrix  :  name of the nifti file output 
%               tICs           :  temporal independent components
%               sICs           :  spatial independent components
%
%  HL 040213
%  Version 1.2

function res=oil_perform_ica(S)

% Set up variables

OSLDIR = getenv('OSLDIR');
% Input data
if isfield(S.results,'ica_concat_path'), ica_concat_path = S.results.ica_concat_path; else, error('No path for ICA input data specified'); end
% ICA settings
if isfield(S,'num_ics'),         num_ics = S.num_ics;              else    num_ics = 20;                        end;
if isfield(S,'last_eig'),        last_eig = S.last_eig;            else    last_eig = 20;                       end;
if isfield(S,'icasso_its'),      icasso_its = S.icasso_its;        else    icasso_its = 0;                      end;
if isfield(S,'nonlinearity'),    nonlinearity = S.nonlinearity;    else    nonlinearity = 'tanh';               end;
if isfield(S,'temp_or_spat'),    temp_or_spat = S.temp_or_spat;    else    temp_or_spat = 'temporal';           end;
if isfield(S,'use_gm_mask'),     use_gm_mask = S.use_gm_mask;      else    use_gm_mask = 0;                     end;
if isfield(S,'gridstep'),        gridstep = S.gridstep;            else    gridstep = 8;                        end;
if isfield(S,'normalise_vox'),   normalise_vox = S.normalise_vox;  else    normalise_vox = 0;                   end; % H Luckhoo and E Hall's recommendation: keep this off. 
if isfield(S,'subj_ind'),        subj_ind = S.subj_ind;            else    subj_ind = [];                       end;
if isfield(S,'initGuess'),       initGuess = S.initGuess;          else    initGuess = 0;                       end;

if use_gm_mask, gm_mask = nii.quickread([OSLDIR '/std_masks/grey_matter/MNI_greymatter_priors_' num2str(gridstep) 'mm.nii.gz'],gridstep);  end

if      icasso_its ~= 0 && strcmp(temp_or_spat,'spatial'),     ica_opt = 'icasso_spatial';
elseif  icasso_its ~= 0 && strcmp(temp_or_spat,'temporal'),    ica_opt = 'icasso_temporal';
elseif  icasso_its == 0 && strcmp(temp_or_spat,'spatial'),     ica_opt = 'fastica_spatial';
else                                                           ica_opt = 'fastica_temporal'; 
end

% Load in concatenated data
% Use weights-normalised data to prevent bias of noise in deep voxels

ica_concat=[];
for i=1:length(ica_concat_path)
    try % mat-file
        [~, ~, fileExt] = fileparts(ica_concat_path{i});
        if strcmpi(fileExt, '.mat'),
            xall=load(ica_concat_path{i});
        else
            xall=load([ica_concat_path{i} '.mat']);
        end%if
        xall=xall.ica_concat;
    catch % assume nifti
        xall=nii.quickread(ica_concat_path{i},gridstep);
    end
    
    if use_gm_mask
        fprintf('\n%s\n','Applying grey matter mask to ICA input');
        maskPrThreshold = 0.3;
        % remove voxels for ica; add them back in later
        xall(gm_mask < maskPrThreshold, :) = [];
    end
    
    ica_concat=[ica_concat; xall];
end


if normalise_vox && ~isempty(subj_ind);
    fprintf('\n%s\n','Normalising all single subject voxels to have unit variance');
    for i=1:length(subj_ind)-1
    ica_concat(:,subj_ind(i):subj_ind(i+1)-1)= normalise(ica_concat(:,subj_ind(i):subj_ind(i+1)-1),2);
    end
end


% Call ICA decomposition

res=S;
switch ica_opt
    case {'icasso_spatial'}
        fprintf('\n%s\n','Performing Spatial ICASSO decomposition');
        if ~initGuess
            fprintf('\n%s\n','Using random initial guess');
            [~,~, res.results.mixing_matrix, res.results.sICs, ~]=icasso(ica_concat',icasso_its,'g',nonlinearity,'lastEig',last_eig,'numOfIC', num_ics,'approach','symm','maxNumIterations',500,'vis','off');
        else
            fprintf('\n%s\n','Using user-provided initial guess');
            [~,~, res.results.mixing_matrix, res.results.sICs, ~]=icasso(ica_concat',icasso_its,'g',nonlinearity,'lastEig',last_eig,'numOfIC', num_ics,'approach','symm','maxNumIterations',500,'vis','off','initGuess',initGuess);
        end
        res.results.mixing_matrix = res.results.mixing_matrix';
    case {'icasso_temporal'}
        fprintf('\n%s\n','Performing Temporal ICASSO decomposition');
        if ~initGuess
            fprintf('\n%s\n','Using random initial guess');
            [~,~, res.results.mixing_matrix, res.results.tICs, ~]=icasso(ica_concat,icasso_its,'g',nonlinearity,'lastEig',last_eig,'numOfIC', num_ics,'approach','symm','maxNumIterations',500,'vis','off');
        else
            fprintf('\n%s\n','Using user-provided initial guess');
            [~,~, res.results.mixing_matrix, res.results.tICs, ~]=icasso(ica_concat,icasso_its,'g',nonlinearity,'lastEig',last_eig,'numOfIC', num_ics,'approach','symm','maxNumIterations',500,'vis','off','initGuess',initGuess);
        end
        res.results.mixing_matrix = res.results.mixing_matrix';
    case {'fastica_spatial'}
        fprintf('\n%s\n','Performing Spatial ICA decomposition');
        if ~initGuess
            fprintf('\n%s\n','Using random initial guess');
            [res.results.sICs,res.results.mixing_matrix,~]=fastica(ica_concat','g',nonlinearity,'lastEig',last_eig, 'numOfIC', num_ics,'approach','symm');
        else
            fprintf('\n%s\n','Using user-provided initial guess');
            [res.results.sICs,res.results.mixing_matrix,~]=fastica(ica_concat','g',nonlinearity,'lastEig',last_eig, 'numOfIC', num_ics,'approach','symm','initGuess',initGuess);
        end
    case {'fastica_temporal'}
        fprintf('\n%s\n','Performing Temporal ICA decomposition');
        if ~initGuess
            fprintf('\n%s\n','Using random initial guess');
            [res.results.tICs,res.results.mixing_matrix,~]=fastica(ica_concat,'g',nonlinearity,'lastEig',last_eig, 'numOfIC', num_ics,'approach','symm');
        else
            fprintf('\n%s\n','Using user-provided initial guess');
            [res.results.tICs,res.results.mixing_matrix,~]=fastica(ica_concat,'g',nonlinearity,'lastEig',last_eig, 'numOfIC', num_ics,'approach','symm','initGuess',initGuess);
        end
    otherwise
end

% Correct the signs of the ICs

switch lower(temp_or_spat)
    case {'temporal'}
        ma=max(res.results.mixing_matrix,[],1);
        mi=min(res.results.mixing_matrix,[],1);
        di=ones(numel(ma),1);
        di(abs(mi)>abs(ma))=-1;
        res.results.mixing_matrix = res.results.mixing_matrix*diag(di);
        res.results.tICs = diag(di)*res.results.tICs;
    case {'spatial'}
        ma=max(res.results.sICs,[],2);
        mi=min(res.results.sICs,[],2);
        di=ones(numel(ma),1);
        di(abs(mi)>abs(ma))=-1;
        res.results.mixing_matrix = res.results.mixing_matrix*diag(di);
        res.results.sICs = diag(di)*res.results.sICs;
        if use_gm_mask,
            tmp = res.results.sICs;
            res.results.sICs = zeros(size(res.results.sICs, 1), size(gm_mask, 1));
            res.results.sICs(:, gm_mask>=maskPrThreshold) = tmp;
        end%if
    otherwise
        error('ICA type not recognised');
end
