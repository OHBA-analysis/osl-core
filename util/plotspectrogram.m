function [varargout] = osl_plotspectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)

[S,F,T,P] = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs);
S = log10(P);


switch nargout,
 case 0,
   imagesc(T,F,S);
   set(gca,'ydir','normal');
   xlabel('time (s)');
   ylabel('frequency (Hz)');
 case 1
   varargout = {S};
 case 2
   varargout = {S,F};
 case 3
   varargout = {S,F,T};
end

end