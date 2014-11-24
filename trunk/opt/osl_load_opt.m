function [opt] = osl_load_opt(optdir)

% load opt structure
%
% [opt] = osl_load_opt(optdir)
% 
% OR:
%
% [opt] = osl_load_opt(opt)
%
%
% MWW 2012

if nargin==0
    [optname, optdir]=uigetfile({'*.mat', 'MAT-files (*.mat)';'*.*', 'All Files'},'Select the opt.mat file');
    opt_fname = fullfile(optdir,optfname);
else
    
    optfname='opt';
    
    if(~isstr(optdir)),
        
        opt=optdir;
        optdir=opt.dirname;        
        
    end;
    
    if(isempty(findstr(optdir, '.opt')))
        optdir=[optdir, '.opt'];
    end;
    
    opt_fname=[optdir '/' optfname '.mat'];
end

tmp=load([opt_fname]);

opt=tmp.opt;
opt.fname=opt_fname;

opt.dirname=optdir;

