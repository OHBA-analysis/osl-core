function [D_tstat, D_cope] = oat_save_spm_stats( Sin )

% [D_tstat, D_cope] = oat_save_spm_stats( Sin )
%
% save out the results of an oat.first_level or oat.group_level as an SPM
% object.
% Requires: Sin.oat,Sin.stats_fname or Sin.stats
%
% Example call:
% 
%   S4=[];
%   S4.oat=oat;
%   S4.stats_fname=oat.group_level.results_fnames;
%   [D_tstat, D_cope]=oat_save_spm_stats(S4); 
%    
% MWW 2013

OSLDIR = getenv('OSLDIR');

if ~isfield(Sin,'suffix')
    Sin.suffix='';
end;

if isfield(Sin,'stats_fname')
    if isstr(Sin.stats_fname')
        Sin.stats=oat_load_results(Sin.oat,Sin.stats_fname);
    else
        Sin.stats=Sin.stats_fname;
    end;
end;

Sin.write_copes=1;

if(Sin.stats.level==1),

    %D=Sin.stats.D_sensor_data;

    D=oat_get_sensordata(Sin.stats);
    
%    if(size(Sin.stats.cope,4)>1)
%        error('Not implemented for time-freq analysis');
%    end;
        
    spm_fname=[Sin.stats.fname Sin.suffix];
        

    if(Sin.write_copes)

        % output nifti file for the contrast of parameter estimates (COPEs)
        % for each contrast

        fnamec=[spm_fname '_cope'  '.dat']; 
        tmp=Sin.stats.cope;       

        Sc=[];
        Sc.D=D;
        Sc.newname=fnamec;
        Sc.newdata=tmp;
        Sc.time=Sin.stats.times;
        Sc.frequencies=Sin.stats.frequencies;
        Sc.modalities=Sin.oat.source_recon.modalities;
        
        for coni=1:length(Sin.stats.contrasts),
            Sc.cond_list{coni}=['cope' num2str(coni)];
        end;
        
        D_cope=osl_change_spm_eeg_data(Sc);
        
        D_cope.save;
    end;

    % output nifti file for the two-tailed tstats for each contrast
    fnamet=[spm_fname '_tstat'  '.dat'];

    tmp=tmp./Sin.stats.stdcope;                        

    Sc=[];
    Sc.D=D;
    Sc.newname=fnamet;
    Sc.newdata=tmp;
    Sc.time=Sin.stats.times;
    Sc.frequencies=Sin.stats.frequencies;
    Sc.modalities=Sin.oat.source_recon.modalities;

    for coni=1:length(Sin.stats.contrasts),
        Sc.cond_list{coni}=['tstat' num2str(coni)];
    end;

    D_tstat=osl_change_spm_eeg_data(Sc);

    D_tstat.save;
     
elseif(Sin.stats.level==2)
     
    % use first subject as the one to show results on
    %D=spm_eeg_load(Sin.oat.source_recon.D_epoched{1}); 
    stats1=oat_load_results(Sin.oat,Sin.oat.first_level.results_fnames{Sin.oat.first_level.sessions_to_do(1)});
    D=oat_get_sensordata(stats1);
    
    for gcon=1:size(Sin.stats.cope,5),
        spm_fname=[Sin.stats.fname Sin.suffix];

        if(Sin.write_copes)

            % output nifti file for the contrast of parameter estimates (COPEs)
            % for each contrast

            tmp=Sin.stats.cope(:,:,:,:,gcon);
            fnamec=[spm_fname '_cope_gc' num2str(gcon) '.dat'];                       

            Sc=[];
            Sc.D=D;
            Sc.newname=fnamec;
            Sc.newdata=tmp;
            Sc.time=Sin.stats.times;
            Sc.frequencies=Sin.stats.frequencies;
            Sc.modalities=Sin.oat.source_recon.modalities;

            for coni=1:size(Sin.stats.cope,3),
                Sc.cond_list{coni}=['cope' num2str(coni)];
            end;

            D_cope{gcon}=osl_change_spm_eeg_data(Sc);
            
            D_cope{gcon}.save;

        end;

        % output nifti file for the two-tailed tstats for each contrast
        fnamet=[spm_fname '_tstat_gc' num2str(gcon) '.dat'];

        tmp=tmp./Sin.stats.stdcope(:,:,:,:,gcon);                        

        Sc=[];
        Sc.D=D;
        Sc.newname=fnamet;
        Sc.newdata=tmp;
        Sc.time=Sin.stats.times;
        Sc.frequencies=Sin.stats.frequencies;
        Sc.modalities=Sin.oat.source_recon.modalities;

        for coni=1:size(Sin.stats.cope,3),
            Sc.cond_list{coni}=['tstat' num2str(coni)];
        end;
            
        D_tstat{gcon}=osl_change_spm_eeg_data(Sc);

        D_tstat{gcon}.save;
    end;
else,
    
    error('level unsupported');
    
end

