close all
clear all

%----------------------------------------------------------------
%%%% PARAMETERS %%%%
%%%%% specify grid size and grid spacing (in m/s):

NX=400; NZ=400; dh=2.0; 		%size of FD grid

%%%% background velocity (in kg/m^3)
vm=3000;

%%%% standard deviation (in percent)
sigma_proz=5.0;

%%%% correlation length (in m):
a=45.0; 

%---------------------------------------------------------------
sigma=sigma_proz*vm/100;


v=fluct(NX,NZ,a,dh,vm,sigma);
 
vm=mean(mean(v));
dev=100*mean(std(v))/vm;

% statistics:
disp(['Maximum Value:',num2str(max(max(v)))]);
disp(['Minimum Value:',num2str(min(min(v)))]);
disp(['Mean Value:',num2str(vm)]);
disp(['Standard deviation (percent):',num2str(dev)]);

x=[1:NX]*dh; 
z=[1:NZ]*dh;
imagesc(x,z,v) 	
xlabel(' X [m]');
ylabel('Y [m]');
colormap(jet)

 
