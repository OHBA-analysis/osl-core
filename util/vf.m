function [tseries,outarg] = vf(Data,varargin)

% alias for ViewFMRI

[tseries,outarg] = feval('ViewFMRI',Data,varargin{:});