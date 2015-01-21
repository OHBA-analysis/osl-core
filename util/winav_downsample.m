function output=winav_downsample(S,dim)

% winav_downsample.m
%
% output = CAE_downsample(S)
%
% Applies windowed averaging to time series designe for estimating 
% the Correlation of Average Envelopes (CAE) described in Brookes 2011.
%
% S must include: 
%   - S.data_in - the array of time courses being downsampled
%   - S.fsample - the sampling frequency of the data in Hz.
%   - S.window - the averaging window length in seconds.
%
% HL - 031011

win=floor(S.window*S.fsample); % CAE window length in samples.

if(win==0), error('win=floor(S.window*S.fsample)=0'); end;

if dim==1
k=1;
from=1;
to=from+win-1;
output=zeros(ceil(size(S.data_in,1)/win)+1,size(S.data_in,2));
while to<=size(S.data_in,1)
    output(k,:,:)=mean(S.data_in(from:to,:,:),1);
    from=to+1;
    to=from+win-1;    
    k=k+1;
end;
output=output(1:k-1,:,:);
end

if dim==2
k=1;
from=1;
to=from+win-1;
output=zeros(size(S.data_in,1),ceil(size(S.data_in,2)/win)+1);
while to<=size(S.data_in,2)
    output(:,k,:)=mean(S.data_in(:,from:to,:),2);
    from=to+1;
    to=from+win-1;
    k=k+1;
end;
output=output(:,1:k-1,:);
end

if dim==3
k=1;
from=1;
to=from+win-1;
dims=size(S.data_in);
dims(dim)=ceil(size(S.data_in,2)/win)+1;
output=zeros(dims);
while to<=size(S.data_in,3)
    output(:,:,k)=mean(S.data_in(:,:,from:to),3);
    from=to+1;
    to=from+win-1;
    k=k+1;
end;
output=output(:,:,1:k-1);
end

