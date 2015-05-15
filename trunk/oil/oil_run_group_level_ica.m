function [ oil ] = oil_run_group_level_ica( oil  )

% [ oil ] = osl_run_group_level_ica ( oil )
%
% takes in an OAT, which needs to be setup by calling oil=osl_setup_oil_for_ica(S), struct 
% and runs group level subjectwise GLM on task positive ICA output.
% 
% This function should normally be called using osl_run_ica(oil);
%
% HL version 1.1 05.02.13


lower_level=oil.ica_first_level;
current_level=oil.ica_group_level;
current_level_results=[];

for sub=1:length(lower_level.results),
   
    % load in single subject stats
    lower_level_results=lower_level.results{sub};
    
    if(sub==1)
        sc=[size(lower_level_results.cope) 1 1];
        cope=zeros(sc(1),length(lower_level.results),sc(2),sc(3),sc(4));
        stdcope=zeros(sc(1),length(lower_level.results),sc(2),sc(3),sc(4));       
    end;

    cope(:,sub,:,:,:)=lower_level_results.cope; % nvox x nsub x ntpts x ncons x nfreqs

    stdcope(:,sub,:,:,:)=lower_level_results.stdcope;
    
end;
%%%%%%%%%%%%%%%%%%%%%%
%% now setup to do GLM stats

num_subjects=size(current_level.group_design_matrix,1);

x=current_level.group_design_matrix;
pinvx=pinv(x);
pinvxtx=pinv(x'*x);

group_contrast=current_level.group_contrast;

if(size(cope,2)~=num_subjects)
    error('mismatched data');
end;

ncons=size(cope,3);
ntpts=1;
nics=oil.ica.num_ics;

% results containers
current_level_results.cope=zeros(nics,ncons,length(group_contrast)); 
current_level_results.stdcope=zeros(nics,ncons,length(group_contrast)); % nics x ncons x ngcons

%% do a 1st level contrast at a time
for c=1:ncons,
    disp(['Computing group statistics']);    
    for vox=1:size(cope(:,:,c),1),    
        dat=cope(vox,:,c)';
        stddat=std(dat);
        dat=dat/stddat;
                            
        [copeout, varcopeout, coapeout, dof]=glm_fast_for_meg(dat,x,pinvxtx,pinvx,group_contrast,0);
        current_level_results.stdcope(vox,c)=sqrt(varcopeout);
        current_level_results.cope(vox,c)=copeout;    
        current_level_results.tstats(vox,c)=copeout./sqrt(varcopeout);
        current_level_results.pvals(vox,c)=t_to_p(current_level_results.tstats(vox,c),num_subjects-1);
    end;
end;


oil.ica_group_level.results=current_level_results;
