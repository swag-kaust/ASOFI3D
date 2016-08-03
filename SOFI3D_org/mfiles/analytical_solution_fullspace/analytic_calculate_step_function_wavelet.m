function [wavelet,wavelet_deriv,wavelet_deriv_2]=calculate_step_function_wavelet(Td,offset,c,npts,delta,Tshift);

% This function calculates the Mueller-Bruestle-function and its first and second time derivative.
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
      wavelet(ii)=0.5-0.75*cos(pi*time_red(ii)/Td)+0.25*cos(pi*time_red(ii)/Td).^3;
   else wavelet(ii)=1;
   end
   if ((time_red(ii) > -delta) && (time_red(ii) < Td+delta))
      wavelet_deriv(ii)=0.75*pi/Td*(sin(pi*time_red(ii)/Td)).^3;
   end
   if ((time_red(ii) > -delta) && (time_red(ii) < Td+delta))
      wavelet_deriv_2(ii)=(9/4)*(pi/Td)^2*(sin(pi*time_red(ii)/Td))^2*cos(pi*time_red(ii)/Td);
   end
end
