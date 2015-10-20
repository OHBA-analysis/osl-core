function C = glean_groupcov(spmfiles)
% Efficiently compute a group covariance from multiple SPM files.
% This does not requiring prior concatenatation of the data.
% See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
%
% C = GLEAN_GROUPCOV(spmfiles)
%
% REQUIRED INPUTS:
%   spmfiles  - A list of SPM12 MEEG objects to compute covariance across
% 
% OUTPUTS:
%   C         - Group covariance matrix
%
% Adam Baker 2015

 
C = 0;
N = 0;
Mu = 0;

for file = spmfiles(:)'
    
    % Compute session-level moments
    D = spm_eeg_load(char(file));
    [c,mu] = glean_cov(D);
    n = sum(~all(badsamples(D,':',':',':')));

    % Update group-level moments:
    C = C + c*(n-1) + (Mu-mu)*(Mu-mu)'*(N*n)/(N+n);
    Mu = (Mu*N + mu*n)/(N + n);
    N = N + n;

end

% Normalise:
C = C ./ (N-1);

end