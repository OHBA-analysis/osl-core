function [ oat ] = oat_consolidate_results( oat )

% [ oat ] = oat_consolidate_results( oat )
%
% Searches the oat directory for files with the appropriate names
% and fills up results_fnames (for each oat stage) for those files 
% that it can find
%
% Note that this overwrites any existing results_fnames

oat=osl_load_oat(oat);

%%%%%%%%%%%%%%%%%%%%%%
%% source recon stage

% check list of SPM MEEG filenames input
if(~isempty(oat.source_recon.D_epoched))
    Ds=oat.source_recon.D_epoched;
elseif(~isempty(oat.source_recon.D_continuous))
    Ds=oat.source_recon.D_continuous;
else
    error('Either oat.source_recon.D_continuous, or oat.source_recon.D_epoched need to be specified');
end
num_sessions=length(Ds);

oat.source_recon.results_fnames=cell(num_sessions,1);
for ii=1:num_sessions
    oat.source_recon.results_fnames{ii}=['session' num2str(ii) '_recon'];
    if ~osl_util.isfile([oat.source_recon.dirname '/' oat.source_recon.results_fnames{ii} '.mat'])
        oat.source_recon.results_fnames{ii}=[];
    end   
end

%%%%%%%%%%%%%%%%%%%%%%
%% first_level stage

oat.first_level.results_fnames=cell(num_sessions,1);
for ii=1:num_sessions
    oat.first_level.results_fnames{ii}=['session' num2str(ii) '_' oat.first_level.name];
    if ~osl_util.isfile([oat.source_recon.dirname '/' oat.first_level.results_fnames{ii} '.mat'])
        oat.first_level.results_fnames{ii}=[];
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%% subject_level stage

num_subjects=length(oat.subject_level.session_index_list);
oat.subject_level.results_fnames=cell(num_subjects,1);
for ii=1:num_subjects
    oat.subject_level.results_fnames{ii}=['subject' num2str(ii) '_' oat.first_level.name '_' oat.subject_level.name];
    if ~osl_util.isfile([oat.source_recon.dirname '/' oat.subject_level.results_fnames{ii} '.mat'])
        oat.subject_level.results_fnames{ii}=[];
    end  
end

%%%%%%%%%%%%%%%%%%%%%%
%% group_level stage

oat.group_level.results_fnames=[oat.first_level.name '_' oat.subject_level.name '_' oat.group_level.name];
if ~osl_util.isfile([oat.source_recon.dirname '/' oat.group_level.results_fnames '.mat'])
    oat.group_level.results_fnames=[];
end
