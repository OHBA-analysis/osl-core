function D = osl_spmfun(fun,S)
% Apply SPM function to MEEG object without creating a new file
% FORMAT D = osl_spmfun(fun,S)
%
% INPUTS:
% fun         - function handle to a valid SPM function
% S           - input structure containing fields necessary for fun
%
% OUTPUTS:
% D           - MEEG object (also written to disk)
%
% Adam Baker 2014

try 
    D = spm_eeg_load(S.D);
catch 
    error('S must contain valid S.D')
end
   
Dnew = feval(fun,S);

if isa(Dnew,'meeg')
    D = move(Dnew,fullfile(D.path,D.fname));
    Dnew.delete;
end

end
