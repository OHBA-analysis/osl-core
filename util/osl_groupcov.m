function C = osl_groupcov(spmfiles)
% Efficiently compute a group covariance from multiple SPM files without
% requiring prior concatenatation of the data.
% See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
%
% Usage:
% C = osl_group_cov(spmfiles)   where spmfiles is a cell array of filenames
%                               for SPM MEEG objects
%
% Adam Baker 2015
%
% TODO: there is currently a small difference when compared to the
% covariance computed over concatenated data - could be precision errors,
% could be a bug...
 
C = 0;
N = 0;
Mu = 0;

for file = spmfiles
    
    % Compute session-level moments
    D = spm_eeg_load(char(file));
    [c,mu] = osl_cov(D);
    n = sum(~all(badsamples(D,':',':',':')));

    % Update group-level moments:
    C = C + c*(n-1) + (Mu-mu)*(Mu-mu)'*(N*n)/(N+n);
    Mu = (Mu*N + mu*n)/(N + n);
    N = N + n;

end

% Normalise:
C = C ./ (N-1);

end



% Check:
% 
% X = {};
% 
% C = 0;
% N = 0;
% Mu = 0;
% 
% for i = 1:10
% 
% X{end+1} = rand(10,750+randi(500));
% 
% n = size(X{end},2);
% c = osl_cov(X{end}) * (n-1);
% mu = mean(X{end},2);
% 
% C = C + c + (Mu-mu)*(Mu-mu)'*(N*n)/(N+n);
% 
% Mu = (Mu*N + mu*n)/(N + n);
% N = N + n;
% 
% end
% 
% C = C ./ (N-1);
% 
% end


