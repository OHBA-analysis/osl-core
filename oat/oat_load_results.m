function res=oat_load_results(oat, fname)
    
% res=oat_load_results.m(oat, fname)
%
% load a matlab struct named results from file named fname

    if(isempty(findstr(oat.source_recon.dirname, '.oat')))
        oat.dirname=[oat.source_recon.dirname, '.oat'];
    end
    
    tmp=load([oat.source_recon.dirname '/' fname]);
    
    res=tmp.oat_stage_results;
    
    if(isfield(res,'source_recon')),
        res.source_recon.dirname=oat.source_recon.dirname;
    end;
    
end

