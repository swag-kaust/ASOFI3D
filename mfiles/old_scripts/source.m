% The function computes amplitude spectrum of a Ricker, Fumue
% or extern source wavelet (central frequency: fs, spectral sampling: df)

function [amp,f]=source(s,fp1,fp2,df,fs)

if ~strcmp(s,'from_file'),
   tsour=1/fs;

   dt=tsour/256;

  k=5;
   tbeg=0;
   tend=k*tsour;

   t=[tbeg:dt:tend];
   n=length(t);
   lsour=tsour/dt;

   if strcmp(s,'fumue'),
      % FUCHS-MUELLER-SIGNAL
      w1=2*pi/tsour;
      ft=sin(w1*(t-tsour))-0.5*sin(2*w1*(t-tsour));
      ft(1:(k/2-0.5)*n/k)=0.0;
      ft((k/2+0.5)*n/k:n)=0.0;
   elseif strcmp(s,'ricker'),
      % RICKER-SIGNAL:
      t0=tsour*1.5;
      T0=tsour*1.5;
      tau=pi*(t-t0)/T0;
      a=4;
      ft=(1-a*tau.*tau).*exp(-2*tau.*tau);
   end

else
   load wavelet.dat
   ft=wavelet(2:size(wavelet,1))';
   dt=1/fs;
end

ft=ft-mean(ft);

% Graphische Darstellung der eingelesenen Zeitreihe
figure;
plot(t,ft)
title(['Eingelesene Zeitreihe']);
xlabel('Zeit [s]');
ylabel('Amplitude');



fnyq=1/(2*dt);
nfft=2^nextpow2(length(ft))

%FFT
y=fftshift(fft(ft,nfft));
amp=abs(y);
amp=amp/max(amp);
f=fnyq*(-nfft/2:nfft/2-1)/(nfft/2);

nn1=(nfft/2+1)+(fp1/fnyq)*nfft/2 ;
nn2=(nfft/2+1)+(fp2/fnyq)*nfft/2;

figure;
plot(f(nn1:nn2),amp(nn1:nn2))
title('Amplitudenspektrum')
xlabel('Frequenz [Hz]')
ylabel('Amplitude');


