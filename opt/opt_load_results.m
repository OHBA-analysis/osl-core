function res=opt_load_results(opt, fname)
    
% res=opt_load_results.m(opt, fname)
%
% load a matlab struct named results from file named fname

    if(isempty(findstr(opt.dirname, '.opt')))
        opt.dirname=[opt.dirname, '.opt'];
    end

    tmp=load([opt.dirname '/' fname]);
    
    res=tmp.opt_results;
    
    if(isfield(res,'logfile')),
        [pth nme ext]=fileparts(res.logfile);
        res.logfile=[opt.dirname '/plots/' nme ext];
    end;
end

