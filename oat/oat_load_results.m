function res=oat_load_results(oat, fname)

% res=oat_load_results.m(oat, fname)
%
% load a matlab struct named results from file named fname

    if(isempty(findstr(oat.source_recon.dirname, '.oat')))
        oat.dirname=[oat.source_recon.dirname, '.oat'];
    end

    % surpress warning if loading older spm8 objects
    warning('off','MATLAB:elementsNowStruc');
    tmp=load([oat.source_recon.dirname '/' fname]);
    warning('on','MATLAB:elementsNowStruc');

    res=tmp.oat_stage_results;

    if(isfield(res,'source_recon')),
        res.source_recon.dirname=oat.source_recon.dirname;
    end;

    % For backwards compatability with osl1.5
    if isfield(res,'mni_coord')
        res.mni_coords = res.mni_coord;
        res = rmfield(res,'mni_coord');
    end

end

