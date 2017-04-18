function Dout=osl_relink_Dinv(Din,oldpath,newpath)
% This little helper function is made to change D.inv paths for structurals
% that are usually hard-coded.
% This is only operating on loaded in D objects. You have to save them by
% Dout.save to have a lasting effect, but keep in mind that you are
% overwriting the original D object then.

% written by RB 2017


Dout=Din;

%oldpath='/home/disk3/ajquinn/Projects/drugface/structurals/M10/';
%newpath='/test/';

% cycles through all invs and changes the path to newpath
for invs=1:length(Din.inv)
    begstr=strfind(Din.inv{invs}.mesh.sMRI,oldpath);
    endstr=length(oldpath)+begstr;
    if ~(isempty(begstr))
        Dout.inv{invs}.mesh.sMRI=prefix(Din.inv{invs}.mesh.sMRI(endstr:end),newpath);
    else
        warning('Path not found')
    end
    
    begstr=strfind(Din.inv{invs}.mesh.def,oldpath);
    endstr=length(oldpath)+begstr;
    if ~(isempty(begstr))
        Dout.inv{invs}.mesh.def=prefix(Din.inv{invs}.mesh.def(endstr:end),newpath);
    else
        warning('Path not found')
    end
    
    begstr=strfind(Din.inv{invs}.mesh.tess_ctx,oldpath);
    endstr=length(oldpath)+begstr;
    if ~(isempty(begstr))
        Dout.inv{invs}.mesh.tess_ctx=prefix(Din.inv{invs}.mesh.tess_ctx(endstr:end),newpath);
    else
        warning('Path not found')
    end
    
    begstr=strfind(Din.inv{invs}.mesh.tess_scalp,oldpath);
    endstr=length(oldpath)+begstr;
    if ~(isempty(begstr))
        Dout.inv{invs}.mesh.tess_scalp=prefix(Din.inv{invs}.mesh.tess_scalp(endstr:end),newpath);
    else
        warning('Path not found')
    end
    
    begstr=strfind(Din.inv{invs}.mesh.tess_oskull,oldpath);
    endstr=length(oldpath)+begstr;
    if ~(isempty(begstr))
        Dout.inv{invs}.mesh.tess_oskull=prefix(Din.inv{invs}.mesh.tess_oskull(endstr:end),newpath);
    else
        warning('Path not found')
    end
    
    begstr=strfind(Din.inv{invs}.mesh.tess_iskull,oldpath);
    endstr=length(oldpath)+begstr;
    if ~(isempty(begstr))
        Dout.inv{invs}.mesh.tess_iskull=prefix(Din.inv{invs}.mesh.tess_iskull(endstr:end),newpath);
    else
        warning('Path not found')
    end
end

