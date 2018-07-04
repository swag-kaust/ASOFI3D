function [wavelet,wavelet_deriv,wavelet_deriv_2]=calculate_sin_hoch_3_wavelet(Td,offset,c,npts,delta,Tshift);

% This function calculates the sin^3-wavelet as a first time derivative of the 
% Mueller-Bruestle-function. Therefore there is a scaling factor included. 
% Also the second time derivative of the Mueller-Bruestle-function is returned by
% this function
%
% Td:     duration of wavelet in s
% offset: offset in m (relevant for calculating seismograms of a pointsource)
% c:      p-wave velocity in m/s (relevant for calculating seismograms of a pointsource) 
% npts:   number of samples
% delta:  sampling interval in s
% Tshift: source onset in s

time=(0:delta:(npts-1)*delta);
time_red=time-offset/c-Tshift;
wavelet=zeros(npts,1);
wavelet_deriv=zeros(npts,1);
wavelet_deriv_2=zeros(npts,1);

for ii=1:length(time)
   if ((time_red(ii) > 0) && (time_red(ii) < Td))
      wavelet(ii)=0.75*pi/Td*(sin(pi*time_red(ii)/Td)).^3;
   end
   if ((time_red(ii) > -delta) && (time_red(ii) < Td+delta))
      wavelet_deriv(ii)=(9/4)*(pi/Td)^2*(sin(pi*time_red(ii)/Td))^2*cos(pi*time_red(ii)/Td);
   end
   if ((time_red(ii) > -delta) && (time_red(ii) < Td+delta))
      wavelet_deriv_2(ii)=(9/4)*(pi/Td)^3*(2*sin(pi*time_red(ii)/Td)*cos(pi*time_red(ii)/Td)^2-sin(pi*time_red(ii)/Td)^3);
   end
end
