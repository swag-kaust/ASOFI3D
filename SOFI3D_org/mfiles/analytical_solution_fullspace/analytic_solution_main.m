%calculated analytic solution of point explosive source in a homogeneous full space


%%%%%%%%%NOTE for comparison with SOFI3D fullspace solution:
% Since the particle velocities are NOT assigned to the GRID points but
% inbetween 2 grid points, the offset of the VX traces in the header is
% slightly incorrect by DH/2, therefore the analytical soultion should be
% calculated with the increment DH/2, in the matlab script
% comparison_analytic_vs_sofi3D.m this offset error is corrected

clear all
close all
clc

%%base filename for output of seismograms
out_file='analyticdata_fullspace'; 

%input parameter
%%%%%%%%%%%%%%%%
vp=3500;                        % p-wave velocity in m/s
vs=2000;                        % s-wave velocity in m/s
rho=2000;                       % density in kg/m^3
npts=1000;                      % number of samples of resulting seismograms
delta=(0.61e-05)*2;             %temporal sampling interval in s of resulting seismograms
dh=0.1/2;                       %spatial sampling interval in m of resulting seismograms

offsets=0.1:dh:(30+dh);         % vector containing the offsets of the seismograms in m
sw_wavelet=2;                   % 1=sin^3-wavelet, else=Ricker
Td=0.002;                       % 1/Td=fc center source frequency
Tshift=0.0*delta;               % source shift in s (as chosen in source.dat)
M0=1;                           % seismic moment in Nm
sw_source_time_function=1;      % 1: wavelet is directly used as source time function (for comparison with  FD solutions)
                                % 2: time integral of wavelet is used as source time function 
DT=delta;                       % sampling interval of SOFI3D modelling 
                                
amp_scale_fact=0.5;             %only for display, factor for scaling the amplitudes


%calculate analytic seismogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (sw_source_time_function==1)
    data=analytic_pointsource_homfullspace(vp,vs,rho,npts,delta,offsets,sw_wavelet,Td,Tshift,M0,sw_source_time_function,DT);
else
    data=analytic_pointsource_homfullspace(vp,vs,rho,npts,delta,offsets,sw_wavelet,Td,Tshift,M0,sw_source_time_function,DT);
end

%export analytic seismogram to SU format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outseis([out_file,'_vx.su'],data,60,80,'vr'); % radial particle velocity
outseis([out_file,'_p.su'],data,60,80,'p'); %pressure seismogram
% outseis('analyticdata_line_ur.su',data,60,80,'ur'); %radial displacement seismogram
%export analytic seismogram to binary format
% eval(['!sustrip < analyticdata_line_p.su > analyticdata_line_p.bin']);


%display analytic seimogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
timescale=delta:delta:delta*npts;

%tracewise normalizing
for ii=1:1:length(offsets)
    data_p(ii,1:npts)=data(ii).trace_p;
    data_vx(ii,1:npts)=data(ii).trace_vr;
    
    data_p(ii,1:npts)=data_p(ii,:)./max(data_p(ii,:)).*amp_scale_fact;
end

% disp(['Maximum VX particle velocity : ',num2str(max(max(data_vx)))]);

hold on
for ii=1:10:length(offsets)
    plot(offsets(ii)+data_p(ii,:),timescale,'LineWidth',2);
end
hold off


title('Analytic Solution');
ylabel('Time in s');
% xlabel('Trace number');
xlabel('Distance in m');
% legend('vert vy','hor vz','Diff');
set(get(gca,'title'),'FontSize',16);
set(get(gca,'title'),'FontWeight','bold');
set(get(gca,'Ylabel'),'FontSize',16);
set(get(gca,'Ylabel'),'FontWeight','bold');
set(get(gca,'Xlabel'),'FontSize',16);
set(get(gca,'Xlabel'),'FontWeight','bold');
set(gca, 'YDir', 'reverse')
set(gca,'FontSize',14);
set(gca,'FontWeight','bold');
set(gca,'Linewidth',1.0);
set(gca,'Box','on');
axis([-1 max(offsets)+1 0 delta*npts]);