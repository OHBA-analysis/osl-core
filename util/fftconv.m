function res=fftconv(x,kern)

% res=fftconv(x,kern)
%
% uses FFT to do convolution
% output is length of x
%
% MWW 2011

if(length(kern)~=length(x))
  error('length(kern)~=length(x)');
end;

res = real(ifft(fft(x).*fft(kern)));
res = res(1:length(x));
