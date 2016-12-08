%display seismogram difference
%compare analytic seismogram and SOFI3D data

%preparation: 
%   1) run the matlab script analytic_solution.m in order to generate analytic data
%   2) run the overnightbuilt test "./benchmark.sh fullspace" in folder ../overnightbuilt in order to generate SOFI3D data
%   3) copy benchmark_fullspace_vx.su and benchmark_fullspace_p.su from folder ../overnightbuilt/su/benchmark_fullspace/ 

%%%%%%%%%NOTE for comparison with SOFI3D fullspace solution:
% Since the particle velocities are NOT assigned to the GRID points but
% inbetween 2 grid points, the offset of the VX traces in the header is
% slightly incorrect by DH/2, therefore the analytical soultion is calculated 
% with the increment DH/2 and later in the script offset index is chosen + 1 
% resutling in [offset + 1/2*DH] 
        
clear all;
close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file names

file_seismo1='./analyticdata_fullspace';
file_seismo2='./benchmark_fullspace';

compare_mode=2; % 1= pressure
                % 2= vx particle velocity 

print2file=1;
filename_out='./SOFI3D_vs_analytic';

%%%%%save picked maximum amplitudes to file

seismo_label1='Analytic solution';
seismo_label2='SOFI3D Data';

choose_offset=[1 5 10 15 20]; %choose trace offset in m for display

%%%%%%normalize trace with respect to 
% 1=global maximum, 
% 2=trace maximum (each input individually), 
% 3=trace max (both input with respect to same maximum),
% else=no
normtrace=3;
amp_scale_fact=2.8;

%%%%%%shift data2 by single samples
shiftdata=0; % 1= automatically determine shift for each trace, then find mean shift value
             % 2= shift data2 by known number of samples (number_shift_samples)
             % else = no
if compare_mode==1
    number_shift_samples=0;
    filename_out=[filename_out,'_p'];
end
if compare_mode==2
    number_shift_samples=1;
    filename_out=[filename_out,'_vx'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loading structured data via su2matlab function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading data');

if compare_mode==1 
    file_seismo1=[file_seismo1,'_p.su'];
    file_seismo2=[file_seismo2,'_p.su'];
end
if compare_mode==2 
    file_seismo1=[file_seismo1,'_vx.su'];
    file_seismo2=[file_seismo2,'_vx.su'];
end

seismo_data1=su2matlab(file_seismo1);
%%%%%%%%%%%%%%%geophone spacing in m
gx1=[seismo_data1.gx]./(-seismo_data1(1).scalco);
offset1=[seismo_data1.offset]./(-seismo_data1(1).scalco);
% gx1=[seismo_data1.gx].*10^(seismo_data1(1).scalco);
% offset1=[seismo_data1.offset].*10^(seismo_data1(1).scalco);

recspac1=(gx1(2)-gx1(1));
%%%%%%%%%%%%%%%sampling rate in s
dt1=seismo_data1(1).d1;
%%%%%sample count
samplespertrace1=seismo_data1.ns;
seismo_data1=[seismo_data1.trace];

seismo_data2=su2matlab(file_seismo2);
gx2=[seismo_data2.gx];
offset2=[seismo_data2.offset]*10^(seismo_data2(1).scalco);
recspac2=(gx2(2)-gx2(1))*10^(seismo_data2(1).scalco);
dt2=seismo_data2(1).d1;
samplespertrace2=seismo_data2.ns;
seismo_data2=[seismo_data2.trace];


%%%%%%cutting upt to maximum time
tmin1 = dt1;
tmin2 = dt2;
tmax = 0.012;


disp('... finished');
disp(' ');

disp('Applying correction, cutting out data, normalizing');
%%%%%%%cutting out time window
if tmax>0
    seismo_data1=seismo_data1(round(tmin1/dt1):round(tmax/dt1),:);
    seismo_data2=seismo_data2(round(tmin2/dt2):round(tmax/dt2),:);
    %%%%%%%creating FD time scales for plotting
    tscale1=(round(tmin1/dt1):1:round(tmax/dt1))*dt1;
    tscale2=(round(tmin2/dt2):1:round(tmax/dt2))*dt2;
else
    %%%%%%%creating FD time scales for plotting
    tscale1=(round(tmin1/dt1):1:samplespertrace1)*dt1;
    tscale2=(round(tmin2/dt2):1:samplespertrace2)*dt2;
end

%determine size of array
xysize1=size(seismo_data1);
xysize2=size(seismo_data2);


%determine shift between data
if shiftdata==1
    for ii=1:1:length(gx1)
        shift_samples(ii)=find(seismo_data1(:,ii)==max(seismo_data1(:,ii)))-find(seismo_data2(:,ii)==max(seismo_data2(:,ii)));
    end
    number_shift_samples=-median(shift_samples);
    shiftdata=2;
end

%shift data2 by known number of samples (number_shift_samples)
if shiftdata==2
    disp(['   shifting data by ',num2str(-number_shift_samples),' sample(s).']);
    seismo_data_tmp=seismo_data2(number_shift_samples+1:end,:);
    clear seismo_data2;
    seismo_data2=seismo_data_tmp;
    clear seismo_data_tmp;
    seismo_data_tmp=seismo_data1(1:end-number_shift_samples,:);
    clear seismo_data1;
    seismo_data1=seismo_data_tmp;
    tscale1=tscale1(1:end-number_shift_samples);
    tscale2=tscale2(1:end-number_shift_samples);
end

disp('... finished');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure1%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting normalized data');


figure(1)

hold on
for ii=1:length(choose_offset)
    
    % determine index corresponding to choosen offset
    if compare_mode==1
        choose_index1=find(offset1==choose_offset(ii));
    elseif compare_mode==2
        %since the particle velocities are NOT assigned to the GRID points but
        %inbetween 2 grid points, the offset in the header is actually
        %wrong, therefore the analytical soultion is calculated with the
        %increment DH/2 and here the offset+1 [offset + 1/2*DH] is chosen
        choose_index1=find(offset1==choose_offset(ii))+1;
    end
            
    choose_index2=find(offset2==choose_offset(ii));
    
    
    % determine factor for normalization
    switch normtrace
        case 1
            % normalizing of the traces with respect to global maximum
            max_amp1=max(max(abs(seismo_data1)));
            max_amp2=max(max(abs(seismo_data2)));
        case 2
            % normalizing of the traces with respect to individual trace maximum
            max_amp1=max(abs(seismo_data1(:,choose_index1)));
            max_amp2=max(abs(seismo_data2(:,choose_index2)));
        case 3
            % normalizing of the traces with respect to trace maximum of input1
            max_amp1=max(abs(seismo_data1(:,choose_index1)));
            max_amp2=max(abs(seismo_data1(:,choose_index1)));
        otherwise
            % no normalizing
            max_amp1=1;
            max_amp2=1;
    end
    
     
    % extracting trace from array amd apply normalizing
    seis_display1=seismo_data1(:,choose_index1)./max_amp1.*amp_scale_fact;
    seis_display2=seismo_data2(:,choose_index2)./max_amp2.*amp_scale_fact;

    % plat traces
    plot(tscale1,offset1(choose_index1)+seis_display1,'-r','LineWidth',2);
    plot(tscale2,offset2(choose_index2)+seis_display2,'--b','LineWidth',2);
    plot(tscale1,offset1(choose_index1)+seis_display1-seis_display2,':g','LineWidth',2);
    
    %calculating root mean square of each trace
    rms(ii,1)=sqrt(1/length(seis_display1)*sum(seis_display1.*seis_display1));
    rms(ii,2)=sqrt(1/length(seis_display2)*sum(seis_display2.*seis_display2));
    rms(ii,3)=sqrt(1/length(seis_display2)*sum((seis_display1-seis_display2).*(seis_display1-seis_display2)));
    disp(['L2 Norm of trace difference of input 1 and input 2 (with respect to input 1, in %) , offset ',num2str(offset1(choose_index1)),' :',num2str(rms(ii,3)/rms(ii,2)*100)]);
    text(0.008, offset1(choose_index1)-1,['RMS Error = ',num2str((rms(ii,3)/rms(ii,2)*100)),' %']);
end
hold off

%%%%%%%%%%%%%%%formatting the figure
ylabel('Offset in m');
xlabel('Time in s');
if compare_mode==1 
    title(['Pressure seismogram']);
elseif compare_mode==2
    title(['v_x particle veloctiy seismogram']);
end

legend(seismo_label1,seismo_label2,'Difference','Location','NorthWest')
set(get(gca,'title'),'FontSize',16);
set(get(gca,'title'),'FontWeight','bold');
set(get(gca,'Ylabel'),'FontSize',16);
set(get(gca,'Ylabel'),'FontWeight','bold');
set(get(gca,'Xlabel'),'FontSize',16);
set(get(gca,'Xlabel'),'FontWeight','bold');

set(gca,'FontSize',14);
set(gca,'FontWeight','bold');
set(gca,'Linewidth',1.0);
set(gca,'Box','on');
axis([0 tmax  -1*amp_scale_fact max(offset1)-5]);


if print2file==1
    savefig(filename_out,'png', '-rgb', '-c0', '-r250');
    eval(['print -depsc  ',[filename_out,'.eps']]);
end
% rms(:,3)./rms(:,1)*100

disp('... finished');
disp('FIN');
