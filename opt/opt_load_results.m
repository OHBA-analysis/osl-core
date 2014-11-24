function res=opt_load_results(opt, fname)
    
% res=opt_load_results.m(opt, fname)
%
% load a matlab struct named results from file named fname

    tmp=load([opt.dirname '/' fname]);
    
    res=tmp.opt_results;
    
end

