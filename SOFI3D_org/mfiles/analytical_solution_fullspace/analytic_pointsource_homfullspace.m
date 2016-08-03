function [data]=calc_pointsource_seismos_homfullspace(vp,vs,rho,npts,delta,offsets,sw_wavelet,Td,Tshift,M0,sw_source_time_function,varargin);

% This function calculates seismograms for an explosive pointsource
% in a homogeneous fullspace.
%
% It is written to compare the solutions of sofi3D (elastic 3D FD forward modelling code) with
% the analytical solution. Therefore one has to keep in mind that the resulting seismograms
% for pressure and velocity/displacement are NOT consistent according to their time axis.
% To obtain an optimal consistence between the pressure seismogram of sofi3D and the
% analytically calculated pressure seismograms one has to shift the analytically determined
% seismograms 1/2*DT (whereas DT is the time sampling interval during FD modelling)
% to smaller times. This has to be done because of the implicit used staggered time grid
% in sofi3D. In this matlab code it is realised by simply shifting the calculated source time
% functions to smaller times.
%
% Input:
%     vp: p-wave velocity in m/s
%     vs: s-wave velocity in m/s
%     rho: density in kg/m^3
%     npts: number of samples of resulting seismograms
%     delta: sampling interval in s of resulting seismograms
%     offsets: vector containing the offsets of the seismograms in m
%     sw_wavelet : switch for wavelet
%                  1: sin^3 impulse
%                  2: Ricker wavelet
%     Td         : duration of source signal in s
%     Tshift     : source onset in s
%     M0         : seismic moment in Nm
%     sw_source_time_function: This switch decides wheather the wavelet given by sw_wavelet is used directly
%                              as the source time function or if the time integral of the wavelet is used as
%                              source time function.
%                              1: wavelet is directly used as source time function
%                              2: time integral of wavelet is used as source time function
%    varargin : This matlab construct is used that you can give an time sampling interval DT used during the
%               FD caluclations. It is only used during this matlab function when sw_source_time_function is
%               set to 2. When no DT is given in this case a warning message is given out that the amplitudes
%               of the calculated seismograms will NOT coincide with the amplitudes of FDMPI seimograms.
%               COMMENT: When the wavelet is not time derivated in the function psource.c in FDMPI but
%                        just added to the elemnts on the main diagonal of the stress tensor (sxx, syy, szz)
%                        without any scaling, then the source time function corresponds to the time integral
%                        of the used wavelet divided by DT (the time sampling interval of the FD calculation).
%
% Output:
%     data:    struct with the following fields
%        .delta    :  sampling interval in s
%        .npts     :  number of samples
%        .offset   :  offset in m
%        .trace_p  :  pressure seismogram
%        .trace_ur :  radial displacement seismogram
%        .trace_vr :  radial velocity seismogram
%


no_of_input_parameters=nargin;

if no_of_input_parameters==12
    DT=varargin{1};
    disp(['DT = ' num2str(DT) 's is used for correcting the']);
    disp('amplitude of the source time function.');
    disp(' ')
elseif (no_of_input_parameters==11 && sw_source_time_function==2)
    disp('The amplitudes of the calculated seismograms will NOT')
    disp('coincide with the amplitudes of SOFI3D seismograms.')
    disp('Furthermore a small time shift will occur between the')
    disp('analytically calculated pressure seismograms and the')
    disp('results of sofi3D.')
    DT=0; % this must be done for the calculation of the source time function for pressure seismograms
end

if ((no_of_input_parameters~=11) && (no_of_input_parameters~=12))
    error('Wrong number of input parameters.')
end


% Calculation of center frequency (more usual when using a Ricker wavelet)
fc=1/Td;

% Definition of time axis
t=(0:delta:(npts-1)*delta);

% Berechnung des Lameschen Parameters lambda
lambda=rho*(vp^2-2*vs^2);
% Berechnung des Kompressionsmoduls
k=(vp^2-(4/3)*vs^2)*rho;



% Calculation of time series
for ii=1:length(offsets)
    data(ii).npts=npts;
    data(ii).delta=delta;
    data(ii).offset=offsets(ii);
    
    % Calculation of the wavelet
    if sw_wavelet==1      % sin^3 impulse
        if sw_source_time_function==1   % directly used as source time function
            
            % source time function for calculating velocity and displacement
            [wavelet,wavelet_deriv,wavelet_deriv_2]=analytic_calculate_sin_hoch_3_wavelet(Td,offsets(ii),vp,npts,delta,Tshift);
            wavelet=(Td/(0.75*pi))*wavelet;
            wavelet_deriv=(Td/(0.75*pi))*wavelet_deriv;
            wavelet_deriv_2=(Td/(0.75*pi))*wavelet_deriv_2;
            
            % source time function for calculating pressure
            [wavelet_p,wavelet_deriv_p,wavelet_deriv_2_p]=analytic_calculate_sin_hoch_3_wavelet(Td,offsets(ii),vp,npts,delta,Tshift-0.5*DT);
            wavelet_p=(Td/(0.75*pi))*wavelet_p;
            wavelet_deriv_p=(Td/(0.75*pi))*wavelet_deriv_p;
            wavelet_deriv_2_p=(Td/(0.75*pi))*wavelet_deriv_2_p;
            
        elseif (sw_source_time_function==2 && no_of_input_parameters==10) % integral used as source time function but no normalization
            
            % source time function for calculating velocity and displacement
            [wavelet,wavelet_deriv,wavelet_deriv_2]=analytic_calculate_sin_hoch_3_wavelet(Td,offsets(ii),vp,npts,delta,Tshift);
            wavelet=(Td/(0.75*pi))*wavelet;
            wavelet_deriv=(Td/(0.75*pi))*wavelet_deriv;
            wavelet_deriv_2=(Td/(0.75*pi))*wavelet_deriv_2;
            
            % source time function for calculating pressure
            [wavelet_p,wavelet_deriv_p,wavelet_deriv_2_p]=analytic_calculate_sin_hoch_3_wavelet(Td,offsets(ii),vp,npts,delta,Tshift-0.5*DT);
            wavelet_p=(Td/(0.75*pi))*wavelet_p;
            wavelet_deriv_p=(Td/(0.75*pi))*wavelet_deriv_p;
            wavelet_deriv_2_p=(Td/(0.75*pi))*wavelet_deriv_2_p;
        else
            
            % source time function for calculating velocity and displacement
            [wavelet,wavelet_deriv,wavelet_deriv_2]=analytic_calculate_step_function_wavelet(Td,offsets(ii),vp,npts,delta,Tshift);
            % Normierung des wavelets
            wavelet=(1/DT)*(Td/(0.75*pi))*wavelet;
            wavelet_deriv=(1/DT)*(Td/(0.75*pi))*wavelet_deriv;
            wavelet_deriv_2=(1/DT)*(Td/(0.75*pi))*wavelet_deriv_2;
            
            % source time function for calculating pressure
            [wavelet_p,wavelet_deriv_p,wavelet_deriv_2_p]=analytic_calculate_step_function_wavelet(Td,offsets(ii),vp,npts,delta,Tshift-0.5*DT);
            % Normierung des wavelets
            wavelet_p=(1/DT)*(Td/(0.75*pi))*wavelet_p;
            wavelet_deriv_p=(1/DT)*(Td/(0.75*pi))*wavelet_deriv_p;
            wavelet_deriv_2_p=(1/DT)*(Td/(0.75*pi))*wavelet_deriv_2_p;
        end
    else                  % Ricker wavelet
        if sw_source_time_function==1  % directly used as source time function
            % source time function for calculating velocity and displacement
            [wavelet,wavelet_deriv,wavelet_deriv_2]=analytic_calculate_Ricker_wavelet(fc,offsets(ii),vp,npts,delta,Tshift);
            
            % source time function for calculating pressure
            [wavelet_p,wavelet_deriv_p,wavelet_deriv_2_p]=analytic_calculate_Ricker_wavelet(fc,offsets(ii),vp,npts,delta,Tshift-0.5*DT);
            
        elseif (sw_source_time_function==2 && no_of_input_parameters==10) % integral used as source time function but no normalization
            % source time function for calculating velocity and displacement
            [wavelet,wavelet_deriv,wavelet_deriv_2]=analytic_calculate_Ricker_wavelet(fc,offsets(ii),vp,npts,delta,Tshift);
            
            % source time function for calculating pressure
            [wavelet_p,wavelet_deriv_p,wavelet_deriv_2_p]=analytic_calculate_Ricker_wavelet(fc,offsets(ii),vp,npts,delta,Tshift-0.5*DT);
        else
            % source time function for calculating velocity and displacement
            [wavelet,wavelet_deriv,wavelet_deriv_2]=analytic_calculate_integral_Ricker_wavelet(fc,offsets(ii),vp,npts,delta,Tshift);
            % Normierung des wavelets
            wavelet=(1/DT)*wavelet;
            wavelet_deriv=(1/DT)*wavelet_deriv;
            wavelet_deriv_2=(1/DT)*wavelet_deriv_2;
            
            % source time function for calculating pressure
            [wavelet_p,wavelet_deriv_p,wavelet_deriv_2_p]=analytic_calculate_integral_Ricker_wavelet(fc,offsets(ii),vp,npts,delta,Tshift-0.5*DT);
            % Normierung des wavelets
            wavelet_p=(1/DT)*wavelet_p;
            wavelet_deriv_p=(1/DT)*wavelet_deriv_p;
            wavelet_deriv_2_p=(1/DT)*wavelet_deriv_2_p;
        end
    end
    
    data(ii).trace_ur=(-M0/(4*pi*vp^2*rho*data(ii).offset^2)).*wavelet-(M0/(4*pi*vp^3*rho*data(ii).offset)).*wavelet_deriv;
    data(ii).trace_vr=(-M0/(4*pi*vp^2*rho*data(ii).offset^2)).*wavelet_deriv-(M0/(4*pi*vp^3*rho*data(ii).offset)).*wavelet_deriv_2;
    %data(ii).trace_p=-k*(2./data(ii).offset).*data(ii).trace_ur-k*(M0/(2*pi*vp^2*rho*data(ii).offset^3))*wavelet-...
    %                                                             k*(M0/(2*pi*vp^3*rho*data(ii).offset^2))*wavelet_deriv-...
    %								k*(M0/(4*pi*vp^4*rho*data(ii).offset))*wavelet_deriv_2;
    data(ii).trace_p=-(k*M0/(4*pi*vp^4*rho*data(ii).offset)).*wavelet_deriv_2_p;
end




