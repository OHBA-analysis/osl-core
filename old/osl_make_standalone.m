function osl_make_standalone(osldir,outdir,toCompile)
% Compile OSL as a standalone executable using the MATLAB compiler
%   http://www.mathworks.com/products/compiler/
%
% This will generate a standalone program, which can be run
% outside MATLAB, and therefore does not use up a MATLAB licence.
%
% Modified from spm_make_standalone.m by gjw 2013

%==========================================================================
%-Start OSL so as to have access to the SPM functions we need
restoredefaultpath;
addpath([osldir '/osl']);
osl_startup(osldir);
%==========================================================================

%==========================================================================
%-Care of startup.m
%==========================================================================
% see http://www.mathworks.com/support/solutions/data/1-QXFMQ.html?1-QXFMQ
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n',...
        fileparts(which('startup')));
end

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spm('dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% create code to insert toolbox config
%-Toolbox autodetection
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d  = dir(tbxdir); d = {d([d.isdir]).name};
dd = regexp(d,'^\.');
%(Beware, regexp returns an array if input cell array is of dim 0 or 1)
if ~iscell(dd), dd = {dd}; end
d  = {'' d{cellfun('isempty',dd)}};
ft = {};
ftc = {};
%-Look for '*_cfg_*.m' or '*_config_*.m' files in these directories
for i=1:length(d)
    d2 = fullfile(tbxdir,d{i});
    di = dir(d2); di = {di(~[di.isdir]).name};
    f2 = regexp(di,'.*_cfg_.*\.m$');
    if ~iscell(f2), f2 = {f2}; end
    fi = {di{~cellfun('isempty',f2)}};
    if ~isempty(fi)
        ft = [ft(:); fi(:)];
    else
        % try *_config_*.m files, if toolbox does not have '*_cfg_*.m' files
        f2 = regexp(di,'.*_config_.*\.m$');
        if ~iscell(f2), f2 = {f2}; end
        fi = {di{~cellfun('isempty',f2)}};
        ftc = [ftc(:); fi(:)];
    end;
end

% get rid of any '._' files
rmv = zeros(size(ft));
for i = 1:length(ft)
    if strcmp('._',ft{i}(1:2))
        rmv(i) = 1;
    end
end
rmv = logical(rmv);
ft(rmv) = [];

if ~isempty(ft)||~isempty(ftc)
    if isempty(ft)
        ftstr = '';
    else
        ft = cellfun(@(cft)strtok(cft,'.'),ft,'UniformOutput',false);
        ftstr  = sprintf('%s ', ft{:});
    end
    if isempty(ftc)
        ftcstr = '';
    else
        ftc = cellfun(@(cftc)strtok(cftc,'.'),ftc,'UniformOutput',false);
        ftcstr = sprintf('cfg_struct2cfg(%s) ', ftc{:});
            end
    fprintf(fid,'values = {%s %s};\n', ftstr, ftcstr);
end

fclose(fid);

%==========================================================================
[compilepath,compileFunction] = fileparts(toCompile);
addpath(compilepath);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
%-Get the path for, and filename of, the function to compile
 
mcc('-m', '-C', '-v', compileFunction,...
    '-d',outdir,...
    '-N','-p',fullfile(matlabroot,'toolbox','signal'),...
    '-p',fullfile(matlabroot,'toolbox','stats'),...
    '-R','-singleCompThread',...
    '-f','./mbuildopts.sh',...
     '-a',osldir);
% 
% ,...
%     '-a',[osldir '/osl'],...
%     '-a',[osldir '/spm8'],...
%     '-a',[osldir '/fmt'],...
%     '-a',[osldir '/myfmt'],...
%     '-a',[osldir '/hmm_box'],...
%     '-a',[osldir '/std_masks'],...
%     '-a',[osldir '/layouts'],...
%     '-a',osldir)

%     '-p',fullfile(matlabroot,'toolbox','stats'),...