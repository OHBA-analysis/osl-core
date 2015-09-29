function [sim,env] = osl_simulate_envelopes(sig,N)
% Simulate multivariate data or envelopes using phase randomisation.
%
% [SIM, ENV] = OSL_SIMULATE_ENVELOPES(SIG,N) returns N realisations of 
% simulated envelopes SIM, which are derived from phase-randomised 
% surrogates of SIG such that any correlation structure in the envelopes of 
% SIG is maintained. Also returns the envelopes of the original data as
% ENV. The input signals SIG should be a [samples x channels] matrix.
% 
% Adam Baker, 2015
%
% See also, phaseran.m by Carlos Gias
%
% Reference:
% Prichard, D., Theiler, J. Generating Surrogate Data for Time Series
% with Several Simultaneously Measured Variables (1994)
% Physical Review Letters, Vol 73, Number 7

% Require odd number of samples for phaseran 
sim = zeros([size(sig),N]);
if mod(size(sim,1),2) == 0
    sim = sim(1:end-1,:,:);
end

% Simulate raw data via phase randomisation
for i = 1:size(sig,2)
    sim(:,i,:) = phaseran(sig(:,i),N);
end

% Compute envelopes of original data and their covariance structure
env = abs(hilbert(sig));   
cholcov = chol(cov(env));

for n = 1:N
    % Compute envelope of simulated data
    sim(:,:,n) = abs(hilbert(sim(:,:,n)));
    % Induce correlation in the envelopes
    sim(:,:,n) = sim(:,:,n) / chol(cov(sim(:,:,n))) * cholcov;
end

% Pad by one sample if necessary
if size(sim,1) ~= size(sig,1)
    sim = [sim; zeros(size(sim(1,:,:)))];
end


end