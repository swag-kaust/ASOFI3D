%creates a sofi3D source file for code-benchmarking
%instead of a real application scenario, mulitple sources are placed
%withing the model domain such as the wavefront excited by each sources
%reaches the neighboring wavenfront at TIME/2
clear all;
close all;
clc;
writefiles=1;%1=yes else=no
srcfilename='../par/sources/sources_codebenchmark.dat';

%input parameters as specified in input file
NX=1024;
NY=1024;
NZ=1024;

DX=1.0;
DY=1.0;
DZ=1.0;

TIME=0.1;

%input parameters as specified in model file
vp=5700;

TD=0.0;
FC=500.0;
AMP=1.0;

%create coordinate vectors

xsrc=round(vp*TIME/(DX*4):vp*TIME/(DX*4):NX);
ysrc=round(vp*TIME/(DY*4):vp*TIME/(DY*4):NY);
zsrc=round(vp*TIME/(DZ*4):vp*TIME/(DZ*4):NZ);

[XSRC,YSRC,ZSRC]=meshgrid(xsrc,ysrc,zsrc);
XSRC=reshape(XSRC,1,[]);
YSRC=reshape(YSRC,1,[]);
ZSRC=reshape(ZSRC,1,[]);

srcfile=zeros(length(xsrc)*length(ysrc)*length(zsrc),6);
tmpvec=ones(length(xsrc)*length(ysrc)*length(zsrc),1);


srcfile(:,1)=XSRC.*DX;
srcfile(:,2)=YSRC.*DY;
srcfile(:,3)=ZSRC.*DZ;
srcfile(:,4)=tmpvec.*TD;
srcfile(:,5)=tmpvec.*FC;
srcfile(:,6)=tmpvec.*AMP;



figure(2);
plot3(srcfile(:,1),srcfile(:,2),srcfile(:,3),'*r');
axis([0 NX*DX 0 NY*DY 0 NZ*DZ])
% set(gca,'YDir','reverse')
set(gca,'DataAspectRatio',[1 1 1]);
box on
xlabel('x in m');
ylabel('y in m');
zlabel('depth z in m');

%write models to disk
disp('...')
if writefiles==1
    file_id=fopen(srcfilename,'w');
    for ii=1:length(srcfile)
        fprintf(file_id,'%2.2f  %2.2f  %2.2f  %2.2f  %2.2f  %2.2f \n',srcfile(ii,1),srcfile(ii,2),srcfile(ii,3),srcfile(ii,4),srcfile(ii,5),srcfile(ii,6));   
    end
    fprintf(file_id,'# \n');   
    fprintf(file_id,'# Definition of (distributed) source positions: \n');
    fprintf(file_id,'# position, time-delay, centre frequency, and amplitude \n');   
    fprintf(file_id,'# \n');   
    fprintf(file_id,'# (comment line is indicated by # or % as first character) \n');   
    fprintf(file_id,'# \n');   
    fprintf(file_id,'# Parameters for each source node (one source node per line): \n');   
    fprintf(file_id,'# \n');   
    fprintf(file_id,'#  XSRC   YSRC   ZSRC   TD   FC   AMP \n');   
    fprintf(file_id,'# \n');
    fprintf(file_id,'# Symbols: \n');
    fprintf(file_id,'# XSRC= x-coordinate of source point [meter] \n');
    fprintf(file_id,'# YSRC= y-coordinate of source point [meter] \n');
    fprintf(file_id,'# ZSRC= z-coordinate of source point [meter] (vertical) \n');
    fprintf(file_id,'# TD= excitation time (time-delay) for source node [s] \n');
    fprintf(file_id,'# FC= centre frequency of source signal [Hz] \n');
    fprintf(file_id,'# AMP= maximum amplitude of source signal \n');
          
    disp(['Total Number of Sources: ',num2str(length(srcfile))]);
    disp(['Source file ',srcfilename,' created and saved to disk']);
    
    fclose(file_id);
end
disp('...')
