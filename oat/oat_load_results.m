function res=osl_load_oat_results(oat, fname)
    
% res=osl_load_oat_results.m(oat, fname)
%
% load a matlab struct named results from file named fname

    tmp=load([oat.source_recon.dirname '/' fname]);
    
    res=tmp.oat_stage_results;
    
    if(isfield(res,'source_recon')),
        res.source_recon.dirname=oat.source_recon.dirname;
    end;
end

