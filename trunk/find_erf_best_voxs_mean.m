function [ res ] = find_erf_best_voxs_new( S )

num_subjects=size(S.cope,2);
tim=S.tim;

copea=S.cope; % num_voxels x num_subject x num_timepoints x num_contrasts

    if(isfield(S,'stdcope'))
        stdcopea=S.stdcope; % num_voxels x num_subject x num_timepoints x num_contrasts
    end;
    
% find best voxel based on z-stats for specified effect contrast
clear y_erf y_best best_voxs best_vox best_score;

con=S.contrast4localisation;

keep=S.num_best_voxs;

% find best vox for each subject
for sub=1:size(copea,2),    

    if(isfield(S,'stdcope'))
        stdcopea=S.stdcope;  
        dat3=permute(copea(:,sub,:,con)./stdcopea(:,sub,:,con),[1 3 2 4]);  
    else
        dat3=permute(copea(:,sub,:,con),[1 3 2 4]);       
    end;

    score=mean(abs(dat3(:,tim>0)),2);
    
    [i,j]=sort(score);

    best_voxs(sub,:)=j(max(1,end-keep+1):end);    

end;
    
res.best_voxs=best_voxs;

end
