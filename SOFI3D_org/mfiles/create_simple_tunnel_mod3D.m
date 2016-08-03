%creates a binary formation model for 3D seismic modeling
close all;
clear all;
clc;

%general parameters
writefiles=1;%1=yes else=no
modname='tunneltest';
plotmodels=1;%1=yes else=no
grid_spacing=0.2;

%dimensions of base model
xmax_base=300;%gridpoints
ymax_base=400; %(vertical)
zmax_base=200;

%tunnel centered in y and z coordinate, parallel to x-axis 
diameter_tunnel=20; %"edge lenght" in grid points -> tunnel is not a tube but a rectangel
%change in material properties, i.e. fault zone
x_start_fault=150; % begin of new material, at fault zone in grid points
y_start_fault=50; % 
z_start_fault=1; % 

%velocity, density base model
vp_base=5700; %m/s
vs_base=3400;
rho_base=2200; %kg/m3

%velocity, density fault zone
vp_fault=4000.0; %m/s
vs_fault=2400.0;
rho_fault=1800.0; %g/cm3

%velocities, density tunnel (assuming vacuum)
vp_tunnel=0.0; %m/s
vs_tunnel=0.00001;
rho_tunnel=1.25; %kg/m3



%----------------VP model
%creating model
mod=(zeros(xmax_base,ymax_base,zmax_base));
%filling in base parmaters
mod(:,:)=vp_base;
%filling in fault parameters (i.e. change in material properties)
mod(x_start_fault:end,y_start_fault:end,z_start_fault:end)=vp_fault;
%filling in tunnel parameters
mod(1:xmax_base,...
    round(ymax_base/2)-round(diameter_tunnel/2):round(ymax_base/2)+round(diameter_tunnel/2),...
    round(zmax_base/2)-round(diameter_tunnel/2):round(zmax_base/2)+round(diameter_tunnel/2))=vp_tunnel;

%plot models
if plotmodels==1
   h1=figure(1);
   subplot(2,2,1);
   pcolor(grid_spacing:grid_spacing:zmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod(xmax_base/2,:,:)));
   shading flat;
   colorbar
   title('vp-velocity Model (y-z-plain)');
   xlabel('z in m (horizontal)');
   ylabel('y in m (vertical)');
   axis equal tight;
   
   subplot(2,2,2);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:zmax_base*grid_spacing,squeeze(mod(:,ymax_base/2,:))');
   shading flat;
   colorbar
   title('vp-velocity Model (x-z-plain)');
   xlabel('x in m');
   ylabel('z in m (horizontal)');
   axis equal tight;
   
   subplot(2,2,3);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod(:,:,zmax_base/2))');
   shading flat;
   colorbar
   title('vp-velocity Model (x-y-plain)');
   xlabel('x in m');
   ylabel('y in m (vertical)');
   axis equal tight;
end

%write models to disk
if writefiles==1
    mod=permute(mod,[2 1 3]);
    fid_mod=fopen(strcat(modname,'.vp'),'w');
    fwrite(fid_mod,mod,'float');
    fclose(fid_mod);
end

clear mod;

%----------------VS model
%creating model
mod=(zeros(xmax_base,ymax_base,zmax_base));
%filling in base parmaters
mod(:,:)=vs_base;
%filling in fault parameters (i.e. change in material properties)
mod(x_start_fault:end,y_start_fault:end,z_start_fault:end)=vs_fault;
%filling in tunnel parameters
mod(1:xmax_base,...
    round(ymax_base/2)-round(diameter_tunnel/2):round(ymax_base/2)+round(diameter_tunnel/2),...
    round(zmax_base/2)-round(diameter_tunnel/2):round(zmax_base/2)+round(diameter_tunnel/2))=vs_tunnel;

%plot models
if plotmodels==1
   h2=figure(2);
   subplot(2,2,1);
   pcolor(grid_spacing:grid_spacing:zmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod(xmax_base/2,:,:)));
   shading flat;
   colorbar
   title('vs-velocity Model (y-z-plain)');
   xlabel('z in m (horizontal)');
   ylabel('y in m (vertical)');
   axis equal tight;
   
   subplot(2,2,2);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:zmax_base*grid_spacing,squeeze(mod(:,ymax_base/2,:))');
   shading flat;
   colorbar
   title('vs-velocity Model (x-z-plain)');
   xlabel('x in m');
   ylabel('z in m (horizontal)');
   axis equal tight;
   
   subplot(2,2,3);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod(:,:,zmax_base/2))');
   shading flat;
   colorbar
   title('vs-velocity Model (x-y-plain)');
   xlabel('x in m');
   ylabel('y in m (vertical)');
   axis equal tight;
end

%write models to disk
if writefiles==1
    mod=permute(mod,[2 1 3]);
    fid_mod=fopen(strcat(modname,'.vs'),'w');
    fwrite(fid_mod,mod,'float');
    fclose(fid_mod);
end

clear mod;

%----------------Rho model
%creating model
mod=(zeros(xmax_base,ymax_base,zmax_base));
%filling in base parmaters
mod(:,:)=rho_base;
%filling in fault parameters (i.e. change in material properties)
mod(x_start_fault:end,y_start_fault:end,z_start_fault:end)=rho_fault;
%filling in tunnel parameters
mod(1:xmax_base,...
    round(ymax_base/2)-round(diameter_tunnel/2):round(ymax_base/2)+round(diameter_tunnel/2),...
    round(zmax_base/2)-round(diameter_tunnel/2):round(zmax_base/2)+round(diameter_tunnel/2))=rho_tunnel;

%plot models
if plotmodels==1
   h3=figure(3);
   subplot(2,2,1);
   pcolor(grid_spacing:grid_spacing:zmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod(xmax_base/2,:,:)));
   shading flat;
   colorbar
   title('Density Model (y-z-plain)');
   xlabel('z in m (horizontal)');
   ylabel('y in m (vertical)');
   axis equal tight;
   
   subplot(2,2,2);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:zmax_base*grid_spacing,squeeze(mod(:,ymax_base/2,:))');
   shading flat;
   colorbar
   title('Density Model (x-z-plain)');
   xlabel('x in m');
   ylabel('z in m (horizontal)');
   axis equal tight;
   
   subplot(2,2,3);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod(:,:,zmax_base/2))');
   shading flat;
   colorbar
   title('Density Model (x-y-plain)');
   xlabel('x in m');
   ylabel('y in m (vertical)');
   axis equal tight;
end

%write models to disk
if writefiles==1
    mod=permute(mod,[2 1 3]);
    fid_mod=fopen(strcat(modname,'.rho'),'w');
    fwrite(fid_mod,mod,'float');
    fclose(fid_mod);
end

clear mod;

