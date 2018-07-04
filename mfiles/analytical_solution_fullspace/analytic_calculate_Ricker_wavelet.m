function [wavelet,wavelet_deriv,wavelet_deriv_2]=calculate_Ricker_wavelet(fc,offset,c,npts,delta,Tshift);

% This function calculates the Ricker wavelet and its first and second time derivative. 
%
% fc:     center frequency in Hz
% offset: offset in m (relevant for calculating seismograms of a pointsource)
% c:      p-wave velocity in m/s (relevant for calculating seismograms of a pointsource) 
% npts:   number of samples
% delta:  sampling interval in s
% Tshift: source onset in s

%time=(0:delta:(npts-1)*delta);
time=(delta:delta:npts*delta);
time_red=time-offset/c-Tshift;
wavelet=zeros(npts,1);
wavelet_deriv=zeros(npts,1);

Td=(1/fc);

tau=pi*(time_red-1.5*Td)/(Td);
wavelet_dummy=(1-2*tau.^2).*exp(-tau.^2);
wavelet_deriv_dummy=(pi/Td)*(4*tau.^3-6*tau).*exp(-tau.^2);
wavelet_deriv_2_dummy=(pi/Td)^2*(-8*tau.^4+24*tau.^2-6).*exp(-tau.^2);
wavelet=wavelet_dummy.';
wavelet_deriv=wavelet_deriv_dummy.';
wavelet_deriv_2=wavelet_deriv_2_dummy.';
clear wavelet_dummy wavelet_deriv_dummy wavelet_deriv_dummy_2 
