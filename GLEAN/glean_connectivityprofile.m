function P = glean_connectivityprofile(hmm)
% Computes the "connectivity profile" P for each HMM state. 
% The connectivity profile is defined for each state as the Euclidean 
% difference between a particular voxel's covariance with all other voxels
% within the state versus outside of the state.
%
% P = GLEAN_CONNECTIVITYPROFILE(hmm)
%
% REQUIRED INPUTS:
%   hmm     - HMM structure inferred by glean_infer_hmm()
% 
% OUTPUTS:
%   P       - [channels x states] connectivity profiles
%
% Adam Baker 2015

Nchannels = size(hmm.state(1).Cov,1);


C = zeros(Nchannels,Nchannels,hmm.K);

for k = 1:hmm.K
    C(:,:,k) = corrcov(hmm.state(k).Cov);
end

P = zeros(Nchannels,hmm.K);

channels = 1:Nchannels;
states = 1:hmm.K;

C(repmat(logical(eye(Nchannels)),[1,1,hmm.K])) = nan;

for channel = channels
    
    for state = states
        
        centroid = mean(C(channel,channels~=channel,states~=state),3);
        cp = C(channel,channels~=channel,state);
        
        P(channel,state) = norm((centroid - cp),2);
        
    end
    
end


% Nchannels = size(hmm.state(1).Cov,1);
% 
% 
% P = zeros(Nchannels,hmm.K);
% 
% for k = 1:hmm.K
%     C = corrcov(hmm.state(k).Cov);
%     P(:,k) = nanmedian(diag(diag(nan(size(C)))) + C);
% end
