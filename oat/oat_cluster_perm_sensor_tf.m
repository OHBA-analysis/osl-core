function [corrp group_tstats dmps] = oat_cluster_perm_sensor_tf(S)

% [corrp group_tstats group_copes] = oat_cluster_perm_sensor_tf(S)
%
% Options:
% S.oat - oat containing group_level oat results to work on
% (oat.group_level.results_fnames)
% S.cluster_stats_thresh - cluster forming thresh
% S.cluster_stats_nperms - bum permutations
% S.modality - modality to run over
% S.first_level_copes_to_do - first level copes to work on
% S.minnbchan - the minimum number of neighbouring channels/combinations
%
% Ouput:
% corrp is cluster corrected P-values
% group_tstats is tstats that were thresholded to form the clusters
%
% MWW 2013

% do 2d or 3d cluster statistics on time-frequency sensor data

OSLDIR = getenv('OSLDIR');

try, masksdir=[OSLDIR '/std_masks' ]; catch, error('OSLDIR not set. Run osl_startup.'); end;

thresh = S.cluster_stats_thresh;
nP = S.cluster_stats_nperms;

try tmp=S.minnbchan; catch, S.minnbchan=0; end;

disp(['Doing cluster perm testing']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load in previously run parametric gstats
gstats=oat_load_results(S.oat,S.oat.group_level.results_fnames);

if ~isfield(gstats,'lower_level_copes'),

    warning('Need lower_level_copes to be stored. Re-running group stage to get them');
    oat=S.oat;
    oat.group_level.store_lower_level_copes=1;
    oat.to_do=[0 0 0 1];
    oat=osl_run_oat(oat);
    S.oat=oat;

    gstats=oat_load_results(S.oat,S.oat.group_level.results_fnames);

end;

S4=[];
S4.oat=S.oat;
S4.stats_fname=S.oat.group_level.results_fnames;
S4.write_copes=1;

[D_tstat, D_cope]=oat_save_spm_stats(S4);

tres=gstats.times(2)-gstats.times(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get modality specific chan layouts

% use first subject as the one to show results on
stats1=oat_load_results(S.oat,S.oat.first_level.results_fnames{S.oat.first_level.sessions_to_do(1)});

switch S.modality,
    case 'MEGMAG',

        layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
        chanind = strmatch('MEGMAG', D_tstat{1}.chantype(stats1.chanind));
        chanindlayout=chanind;
        chanindoat = strmatch('MEGMAG', D_tstat{1}.chantype(stats1.chanind));

    case 'MEGPLANAR',

        layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
        chanind = strmatch('MEGCMB', D_tstat{1}.chantype(stats1.chanind));
        chanindlayout=chanind;
        chanindoat = strmatch('MEGCMB', D_tstat{1}.chantype(stats1.chanind));

    case 'MEG',

        layout      = [OSLDIR '/layouts/CTF275.lay'];
        chanind = strmatch('MEG', D_tstat{1}.chantype(stats1.chanind));
        chanindlayout=chanind;
        chanindoat = strmatch('MEG', D_tstat{1}.chantype(stats1.chanind));

    case 'EEG',

        B=sensors(D_tstat{1},'EEG');
        mat=spm2fieldtrip(D_tstat{1});
        cfgx=[];
        cfgx.elec.pnt=B.elecpos;
        cfgx.elec.label=B.label;
        lay = ft_prepare_layout(cfgx,mat);

        layout =lay;
        chanind = strmatch('EEG', D_tstat{1}.chantype(stats1.chanind));
        chanindlayout=chanind;
        chanindoat = strmatch('EEG', D_tstat{1}.chantype(stats1.chanind));

end;
clear chanind;

data.label=D_tstat{1}.chanlabels(chanindlayout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make sensor neighbour structure and check it
cfg_nei = [];
cfg_nei.method = 'triangulation';
% cfg.neighbourdist = 0.4;
cfg_nei.layout = layout;

neighbours = ft_prepare_neighbours(cfg_nei);

cfg_nplot = [];
cfg_nplot.neighbours = neighbours;
cfg_nplot.layout  = layout;
ft_neighbourplot(cfg_nplot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over first level contrasts
current_level=S.oat.group_level;

for coni=1:length(S.first_level_copes_to_do),

    con=S.first_level_copes_to_do(coni);

    cope_smooth_lower_level=gstats.lower_level_copes{con};

    c=permute(cope_smooth_lower_level(chanindoat,:,:,:),[2 1 4 3]);

    nS = size(c,1); %number of 'subjects' (recording sites)
    nR = size(c,2); %number of sensors
    nF = size(c,3); %number of frequencies
    nT = size(c,4); %number of timebins

    cr = reshape(c,nS,nR*nF*nT);

    %build up a null distribution of the maximum cluster size for each
    %permuatation i and each regressor r:

    % setup gauss for temporal smoothing
    if(S.group_varcope_time_smooth_std>0)
        fft_gauss=fftshift(gauss(S.group_varcope_time_smooth_std/tres,1,nT)');
    end;

    dm=S.oat.group_level.group_design_matrix';

    for gconi=1:length(S.group_level_copes_to_do),

        gcon=S.group_level_copes_to_do(gconi);

        % establish type of perms needed
        gcontrast=S.oat.group_level.group_contrast{gcon};
        do_signflip=0;
        do_rowpermute=0;
        str=['Group contrast:' mat2str(gcontrast)];

        c=S.oat.group_level.group_contrast{gcon};
        Q=pinv(dm'*dm);
        F=pinv(c'*Q*c);
        Xeff= dm*Q*c*F; % see effective regressor in appendix of Smith et al. Meaningful design and contrast estimability in FMRI. Neuroimage (2007) vol. 34 (1) pp. 127-36

        if(range(sign(Xeff))==0)
           do_signflip=1; % constant regressor
        else
           do_rowpermute=1;
        end;

        for gg=1:length( gcontrast ),
            if(gcontrast(gg)~=0),
                if(range(dm(:,gg))==0)
                   do_signflip=1; % constant regressor
                else
                   do_rowpermute=1;
                end;
            end;
        end;
        if(do_signflip), str=[str ', doing sign flipping']; end;
        if(do_rowpermute), str=[str ', doing row permuting']; end;

        disp(str);

        ft_progress('init', 'etf');

        % loop over permutations
        for i = 1:nP,
            ft_progress(i/nP);

            % get permuted design matrix
            dmp=permute_dm(dm,gcontrast,do_signflip,do_rowpermute);
            dmps{i}=dmp;

            [cg,vg,tg] = ols(cr,dmp,gcontrast');

            % variance smooth over time:
            if(S.group_varcope_time_smooth_std>0)
                vg=reshape(vg,1,nR,nF,nT);
                for vox=1:size(2,1),
                    for f=1:size(vg,3),
                        dat=permute(vg(1,vox,f,:),[4, 1, 2, 3]);
                        dat2 = fftconv(dat,fft_gauss);
                        vg(1,vox,f,:)=dat2;
                    end;
                end;

                cg=reshape(cg,1,nR,nF,nT);
                tg=cg./sqrt(vg);

            else
                tg = reshape(tg,1,nR,nF,nT);
            end;

            neighbours_dim=2;
    %        [imlabel LL] = spm_bwlabel(double(permute(tg(1,:,:,:),[2 3 4 1])>thresh),26);
            [imlabel LL] = label_clusters(double(permute(tg(1,:,:,:),[2 3 4 1])>thresh),S.minnbchan,neighbours,neighbours_dim,data.label);

            tmp = unique(imlabel);
            tmp(tmp==0) = [];
            nL = 0;
            if ~isempty(tmp)
              for k = tmp' %loop over clusters
            nL(k) = sum(sum(imlabel==k));
              end
            end
            nulldist(i)=max(nL);

        end;

        ft_progress('close');

        %run a one sample t-test on data and compare cluster size to
        %permutation data
        [cg,vg,tg] = ols(cr,dm,gcontrast');

        % variance smooth over time:
        if(S.group_varcope_time_smooth_std>0)
            vg=reshape(vg,1,nR,nF,nT);
            for vox=1:size(vg,2),
                for f=1:size(vg,3),
                    dat=permute(vg(1,vox,f,:),[4, 1, 2, 3]);
                    dat2 = fftconv(dat,fft_gauss);
                    vg(1,vox,f,:)=dat2;
                    %figure;plot(dat);ho;plot(dat2,'r');
                end;
            end;

            cg=reshape(cg,1,nR,nF,nT);
            tg=cg./sqrt(vg);

        else
            tg = reshape(tg,1,nR,nF,nT);
        end;

        neighbours_dim=2;
    %    [imlabel LL] = spm_bwlabel(double(squeeze(tg(1,:,:,:))>thresh),26);
        [imlabel LL] = label_clusters(double(permute(tg(1,:,:,:),[2 3 4 1])>thresh),S.minnbchan,neighbours,neighbours_dim,data.label);

        tmp = unique(imlabel);
        tmp(tmp==0) = [];
        cpimlabel = zeros(size(imlabel));
        if ~isempty(tmp)
        for k= tmp' %loop over clusters
          nL = sum(sum(imlabel==k));
          cp = mean(nL>nulldist);
          cpimlabel(imlabel==k)=cp;
        end
        end

        corrp(chanindlayout,:,con,:,gcon) = cpimlabel;
        group_tstats(chanindlayout,:,con,:,gcon) = permute(tg(1,:,:,:),[2 4 3 1]);
    %    group_copes(chanindoat,:,con,:,gcon) = permute(cg(1,:,:,:),[2 4 3 1]);

    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [imlabel LL] = label_clusters(bw, minnbchan, neighbours, neighbours_dim, labels)

% Input:
% bw:         : (Binary) image to perform labelling on. Can
%               be 2 or 3D.
% minnbchan   :  the minimum number of neighbouring channels/combinations
%
% Output:
% imlabel     : Connected component image, i.e. image where
%               each non-zero voxel in bw will have a value
%               corresponding to its label.
% LL          : Number of connected components in L.

cfg=[];
cfg.neighbours=neighbours;
cfg.channel=labels;
cfg.avgoverchan='no';
channeighbstructmat = makechanneighbstructmat(cfg);

dims=size(bw);
[imlabel LL] = findcluster(reshape(bw, [dims(1) prod(dims(2:end))]),channeighbstructmat,minnbchan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channeighbstructmat = makechanneighbstructmat(cfg);

% MAKECHANNEIGHBSTRUCTMAT makes the makes the matrix containing the channel
% neighbourhood structure.

% because clusterstat has no access to the actual data (containing data.label), this workaround is required
% cfg.neighbours is cleared here because it is not done where avgoverchan is effectuated (it should actually be changed there)
% Adapted from inside fieldtrip func clusterstat.m - MWW 2013

if strcmp(cfg.avgoverchan, 'no')
  nchan=length(cfg.channel);
elseif strcmp(cfg.avgoverchan, 'yes')
  nchan = 1;
  cfg.neighbours = [];
end
channeighbstructmat = false(nchan,nchan);
for chan=1:length(cfg.neighbours)
  [seld] = match_str(cfg.channel, cfg.neighbours(chan).label);
  [seln] = match_str(cfg.channel, cfg.neighbours(chan).neighblabel);
  if isempty(seld)
    % this channel was not present in the data
    continue;
  else
    % add the neighbours of this channel to the matrix
    channeighbstructmat(seld, seln) = true;
  end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dmp=permute_dm(dm,gcontrast,do_signflip,do_rowpermute)

dmp=dm;
if(do_signflip),
    tmp=repmat(((rand(size(dm,1),1)>0.5)-0.5)*2,1,size(dm,2));
    dmp = dm.*tmp;
else
    dmp = dm;
end;

if(do_rowpermute),
    loop=1;
    while loop
        dmp2=dmp(randperm(size(dmp,1)),:);
        sum(abs(gcontrast'*dmp2'-gcontrast'*dmp'))

        if sum(abs(gcontrast'*dmp2'-gcontrast'*dmp'))~=0
            loop=0;
        else
            'same'
        end;
    end;
    dmp=dmp2;
end;




