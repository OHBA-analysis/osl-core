function bf_save(BF, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_save.m 104 2014-01-15 17:45:40Z litvak.vladimir@gmail.com $
%
% This is a patch that allows us saving in -v7, otherwise incredibly bad performance 
% System: MBP '15, Matlab R2016a, OSx 10.11.6 
% RB 2017

if nargin == 1 && exist(fullfile(pwd, 'BF.mat'), 'file')
    save('BF.mat', '-struct', 'BF', '-append');
else
    %    save('BF.mat', '-struct', 'BF', '-v7.3');
    
    % appending to -v7.3 slows down performance incredibly, so modified here
    % todo: is that sufficient storage capacity for us?
    save('BF.mat', '-struct', 'BF', '-v7');
    
end