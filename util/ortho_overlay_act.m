function o = ortho_overlay_act( S )
                
% ortho_overlay_act( S )
%
% S.fname is name of nii file in MNI space
% S.vol_index is vol index (counting from 1) - if not set will assume fname
% is 3D
% S.mni_coord in mm
% S.gridstep is desired spatial resolution in mm
% S.title
% S.percrange (sets S.range using percentiles) 
% or S.range (sets range directly)

try, tmp=S.gridstep; catch, S.gridstep=2; end;

fname=S.fname;

% get current gridstep
[ mni_res ] = nii.get_spatial_res( fname );
mni_res=mni_res(1);

% resample volume
if mni_res~=S.gridstep,
    [pth name ext]=fileparts(fname);
    [~,name] = fileparts(name); % This is an bit of a hack to handle .nii.gz
    new_fname = fullfile(pth,sprintf('%s_%dmm.nii.gz',name,S.gridstep));
    
    if ~isfield(S,'interp')
        S.interp = 'cubic';
    end

    tmp = nii.resample(fname, new_fname, S.gridstep, 'interptype',S.interp);
    
    fname=new_fname;
    mni_res=S.gridstep;
end;

% find index for mni coord
map=nii.load(fname);
map=mean(map,4); 
x1 = abs(map(map~=0));

if isempty(x1)
    low = [];
    high = [];
else
    if isfield(S, 'percrange')
        try
            low=percentile((x1),S.percrange(1));
            high=percentile((x1),S.percrange(2));
        catch
            low=min(x1);
            high=max(x1);
        end
    else
        low=S.range(1);
        high=S.range(2);
    end
    
    if(low==high) % Just use auto limits if this occurs
        low=[];
        high=[];
    end
end

o = osleyes(fname,'show_controls',false,'show_crosshair',false,'title',S.title,'current_point',S.mni_coord,'clim',{[],[low high]});

if isfield(S,'vol_index')
    o.current_vols(2) = S.vol_index;
end

if isfield(S,'title')
    o.title = S.title;
end


end
