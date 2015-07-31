function P = glean_connectivityprofile(hmm)
% Computes the connectivity profile P for each HMM state
%
% P = glean_connectivityprofile(hmm)

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
