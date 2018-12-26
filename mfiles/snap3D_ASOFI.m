%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---script for the visualization of snapshots gained from the ASOFI simulation
%---most parameters are as specified in ASOFI parameter-file, e.g. sofi3D.json
%---Please note : y denotes the vertical axis!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all; 
clearvars; clc;

addpath('./utils');

% User-defined parameters.
% Directory name with the simulation input and output, relative to this script.
plot_opts.par_folder = '../par';
% % Path to configuration file, relative to par_folder.
plot_opts.config_file='./in_and_out/sofi3D.json';
plot_opts.file_out = [plot_opts.par_folder, '/figures/'];
plot_opts.file_ext = '.bin.div';

for phi2=0:15:90
    plot_opts.phi2 = phi2;
    snap3D_ASOFI_fun(plot_opts);
end
 
function snap3D_ASOFI_fun(plot_opts)
phi2 = plot_opts.phi2;
par_folder = plot_opts.par_folder;
config_file = plot_opts.config_file;


%% Read from json to opts.
opts = read_asofi3D_json([par_folder, '/', config_file]);

%% merge snapshots if they were not merged before

oldpwd = pwd;
cd(par_folder)
snap_name_full = [opts.SNAP_FILE,'.bin.div'];
dir_Full = dir(snap_name_full);
dir_000 = dir([opts.SNAP_FILE,'.bin.div.0.0.0']);

if ~exist(snap_name_full,'file')
    system(['../bin/snapmerge ' config_file]);
elseif dir_Full.datenum < dir_000.datenum
    system(['../bin/snapmerge ' config_file]);
end

cd(oldpwd)

%% prepare colormap
create_colormaps;

%% read parameters from json
nx = str2num(opts.NX);
ny = str2num(opts.NY);
nz = str2num(opts.NZ);

outx = str2num(opts.IDX);
outy = str2num(opts.IDY);
outz = str2num(opts.IDZ);

dh = str2num(opts.DX); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---input, output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input file1 (snapshot file1)
file_inp1 = [par_folder,'/snap/test.bin.div'];

% Model file (for single display or contour plot ontop of snapshot)
file_mod = [par_folder,'/model/test.SOFI3D.rho'];

% Output file
% switch for saving snapshots to picture file 1=yes (jpg) 2= yes (png) other=no
filesave=1;
% base name of picture file output, will be expanded by extension jpg/png
file_out=plot_opts.file_out;



% title strings for each sub-figure
title_inp1='\nabla \cdot u';

title_mod='Density model';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---variety of switches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch for contour of model overlaying the model or snapshot file
% 1=yes other=no
cont_switch=1;
% number of contours for the contour-plot
numbOFcont=8;

% Choose model or snapshot plot (model-->1; snapshot-->2)
type_switch=2;

% Choose slice geometry and postion (for 2-D plots)
slice_switch=1; % horizontal(zx)=1; vertical (yx)=2; vertical (yz)=3;
% slice definition, where to slice trough the 3-D space
nx=nx/outx;ny=ny/outy;nz=nz/outz;
zslice=nz/2; % for xy plane
yslice=ny/2; % for xz plane
xslice=nx/2; % for yz plane in grid points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Snapshot definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time increment for snapshots:
TSNAP1=str2num(opts.TSNAP1);
TSNAPINC=str2num(opts.TSNAPINC);
% firts and last snapshot that is considered for displayin
firstframe=2;
lastframe=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---3D definitions: defines two rotating planes (xz, yz plane)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi1=0; % a horizontal plane (x-z plane) is rotated by phi1 with respect to the rotpoint
%phi2=90; % a horizontal plane (x-z plane) is rotated by phi2 with respect to the rotpoint
% rotaxis and rotpoint refers to the rotation of the 2D-slices within the 3D volume
% direction rotation axes [0,1,0] rotation of plane around vertical axis
% [1,0,0] rotation of plane around x-axis
rotaxis=[0,1,0];
% defines point for rotation of the 2D slices [x y z] in meter
%  values are defined as difference to the center of the model/snaphot
% in case of rotpoint=[0,0,0], this point is excatly the center of the model/snaphot
rotpoint=[0 0 0];
% defines angles of perspective for the 3-D view
% i.e. it rotates the plot to favorable perspective
viewpoint=[1,4,1];

file_out = [file_out, num2str(phi2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---axis limits for 2D and 3D visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% colorbar boundaries for cropping the snapshot value range
% only used if type_switch=2
auto_scaling=1; % 1= automatic determination of boundaries, 2= constant values caxis_value , 3= no scaling
caxis_value_1=5e-11;


% delay between the display of each snaphot in s
pause4display=0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---end of input parameter definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---creation of model vectors and loading file data----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot range and increment
xp1=dh; xp2=nx*dh; yp1=dh; yp2=ny*dh; zp1=dh; zp2=nz*dh;

% Computing range for axis and subscript range for the movie
x=xp1:dh*outx:xp2*outx;
y=yp1:dh*outy:yp2*outy; % vertical axis
z=zp1:dh*outz:zp2*outz;

% if model is plotted, than there is only one snapshot
% loop over # of snapshots will terminate after one iteration
if type_switch==1
    firstframe=1
    lastframe=1;
end

% load model file ; Format ieee-be for Sun or Workstation, -le for PC
if (cont_switch==1) | (type_switch==1)
    
    % opening file and reading
    disp(['Loading file ' file_mod]);
    [fid_mod, err_msg] = fopen(file_mod, 'r', 'ieee-le');
    if fid_mod == -1
        disp(['ERROR: Cannot open file ' file_mod]);
        disp(['Reason: ' err_msg]);
        return
    end
    mod_data=fread(fid_mod,'float');
    mod_data=reshape(mod_data,ny,nx,nz);
    mod_data=permute(mod_data,[2,3,1]);
    fclose(fid_mod);
    
    if (slice_switch==1) % z-x-plane (horizontal)
        mod_data=squeeze(mod_data(:,:,yslice));
        mod_data=permute(mod_data,[2,1]);
    end
    if (slice_switch==2) % y-x-plane (vertical)
        mod_data=squeeze(mod_data(:,zslice,:));
        mod_data=permute(mod_data,[2,1]);
    end
    if (slice_switch==3) % y-z-plane (vertical)
        mod_data=squeeze(mod_data(xslice,:,:));
        mod_data=permute(mod_data,[2,1]);
    end
    
end

% open snapshot data of 1st input file
disp(['Loading snap shot file ' file_inp1]);
fid_file1=fopen(file_inp1,'r','ieee-le');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---3D display of snapshot data (single snapshot only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if type_switch==1
        % loading model data
        fid=fopen(file_mod,'r','ieee-le');
        lastframe=1;
        
        if cont_switch==1
            fid_mod=fopen(file_mod,'r','ieee-le');
            mod_data=fread(fid_mod,(nx*ny*nz),'float');
            mod_data=reshape(mod_data,ny,nx,nz);
            mod_data=permute(mod_data,[3,2,1]);
            fclose(fid_mod);
        end
    end
    if type_switch==2
        % loading snapshot data of input1
        fid=fopen(file_inp1,'r','ieee-le');
        colormap(load('./srgb.map'));
        
        % adding contour of model to snapshot data
        if cont_switch==1
            fid_mod=fopen(file_mod,'r','ieee-le');
            mod_data=fread(fid_mod,(nx*ny*nz),'float');
            mod_data=reshape(mod_data,ny,nx,nz);
            mod_data=permute(mod_data,[3,2,1]);
            fclose(fid_mod);
        end
    end
    
    % loop over number of snaphsots
    for i=firstframe:lastframe
        
        % calculating time of snapshot
        tsnap=(i-1)*TSNAPINC+TSNAP1;
        disp(['Loading snapshot no ',int2str(i),' at time=',num2str(tsnap),' s.']);
            
        % loading data (model or snapshot)
        % calculate offset in order to jump to specific snapshot within file
        offset=4*nx*ny*nz*(i-1);
        fseek(fid,offset,-1);
        file1_data=fread(fid,(nx*ny*nz),'float');
        file1_data=reshape(file1_data,ny,nx,nz);
        file1_data=permute(file1_data,[3,2,1]);
        
        D = merge_snapshots(par_folder, plot_opts.file_ext);
        file1_data = D(:,:,:,i);
        % creating a grid
        [X,Z,Y]=meshgrid(x,z,y);
        
        %surface height of horizontal y-x plane
        xyplane=zeros(ny,nx);
        xyplane(:,:)=(nz*dh*outy)/2+rotpoint(2); 
                      
        % creating a vertical slice plane in y-x plane
        hslice = surf(x,y,xyplane);
        % rotate slice plane with, with respect to a point from which rotation is defined
        % in case of rotpoint=[0,0,0], this point is the center of the model/snaphot
        rotate(hslice,rotaxis,phi1,[mean(x)-rotpoint(1),mean(y)-rotpoint(2),mean(mean(xyplane))-rotpoint(3)]);
        % get boundaries of rotated plane
        xd = get(hslice,'XData');
        yd = get(hslice,'YData');
        zd = get(hslice,'ZData');
        % remove plane, only limits are further used
        delete(hslice);
        
        % creating a horizontal slice plane in z-x plane
        hslice2 = surf(x,y,xyplane);
        % rotate slice plane with, with respect to a point from which rotation is defined
        % in case of rotpoint=[0,0,0], this point is the center of the model/snaphot
        rotate(hslice2,rotaxis,phi2,[mean(x)-rotpoint(1),mean(y)-rotpoint(2),mean(mean(xyplane))-rotpoint(3)]);
        % get boundaries of rotated plane
        xd2 = get(hslice2,'XData');
        yd2 = get(hslice2,'YData');
        zd2 = get(hslice2,'ZData');
        % remove plane, only limits are further used
        delete(hslice2);
        
        % display sliced and rotated plane
        h1=figure(1);
        % determination of screen size
        scrsz = get(0,'ScreenSize');
        % determination size of window for plotting
        %set(h1,'Position',[1 scrsz(4)*2/3 scrsz(3)*1/4 scrsz(4)*2/3]);
        axis equal
        % strict vertical slice, no rotation applied
%                 h = slice(Y,X,Z,file1_data,[],yslice*dh*outy-1,[]); 
        % rotated vertical slice
        % !!! note that taking sclices from a homogeneous model is - for some
        % reason not working, please use 2-D visualization instead !!!
        h = slice(X,Z,Y,file1_data,xd,zd,yd);
        set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha',0.5);
        
        hold on
        % strict horizontal, no rotation applied
        %         h2 = slice(Y,X,Z,file1_data,[],[],zslice*dh*outz); % this is strict vertical, no rotation applied
        % rotated horizontal slice
        h2 = slice(X,Z,Y,file1_data,xd2,zd2,yd2);
        set(h2,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha',0.5);
        axis equal
        
        
        if cont_switch==1
            % vertical contour slice
            h=contourslice(X,Z,Y,mod_data,xd,zd,yd,numbOFcont);
            set(h,'FaceColor','black','EdgeColor','black');
            % horizontal contour slice
            h2=contourslice(X,Z,Y,mod_data,xd2,zd2,yd2,numbOFcont);
            set(h2,'FaceColor','black','EdgeColor','black');
        end
        hold off
        
        % formating figure
        %daspect([1,1,1]);
        axis tight
        box on
        
        % set viewpoint within 3-D space
        % (rotate 3-D figure plot to favorable perspective)
        view(viewpoint);
        camproj perspective
        % lightangle(-45,45);
        
        set(gcf,'Renderer','zbuffer');
        %colorbar
        xlabel('x in m');
        ylabel('y in m');
        zlabel('Depth z in m');
        
        % adding title string to plot
        if type_switch==2
            title(title_inp1);
        end
        if type_switch==1
            title(title_mod);
        end
        
        set(get(gca,'title'),'FontSize',13);
        set(get(gca,'title'),'FontWeight','bold');
        set(get(gca,'Ylabel'),'FontSize',12);
        set(get(gca,'Ylabel'),'FontWeight','bold');
        set(get(gca,'Xlabel'),'FontSize',12);
        set(get(gca,'Xlabel'),'FontWeight','bold');
        set(get(gca,'Zlabel'),'FontSize',12);
        set(get(gca,'Zlabel'),'FontWeight','bold');
        set(gca, 'ZDir', 'reverse')
        set(gca, 'YDir','reverse')
        set(gca,'FontSize',12);
        set(gca,'FontWeight','bold');
        set(gca,'Linewidth',1.0);
        set(gca,'Box','on');
        
        % determing maximum amplitude of plot
        file1max=max(max(max(abs(file1_data))));
        % display maximum amplitude to command window
        disp(['  Maximum amplitude of ',title_inp1,'-snapshot: ', num2str(file1max)]);
        
        % cropping model/snapshot dimensions
        ylim([0 nz*dh*outz+0.1])
        xlim([0 nx*dh*outx+0.1])
        zlim([0 ny*dh*outy+1.1])
        
        % for snapshot input files
        if type_switch==2
            
            % limiting the colorbar to specified range
            switch auto_scaling
                case 1
                    caxis([-file1max/5 file1max/5]);
                case 2
                    caxis([-caxis_value_1 caxis_value_1]);
                otherwise
            end
            
            % adding text string for timestep of snapshot
            set(text(3500,-1000,['T= ',sprintf('%2.4f',tsnap), 's']),'FontSize',12,'FontWeight','bold','Color','black');
            
            pause(pause4display);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---Saving the snapshot to file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (filesave~=0)
                % generating filename string
                if i<10
                    imagefile=[file_out,'3D_00',int2str(i)];
                else if i<100
                        imagefile=[file_out,'3D_0',int2str(i)];
                    else
                        imagefile=[file_out,'3D_',int2str(i)];
                    end
                end
                
                % output as jpg graphic (or eps) via print command
                if filesave==1
                    eval(['print -djpeg100  ' [imagefile,'.jpg']]);
                    %             eval(['print -depsc  '[imagefile2,'.eps']]);
                end
                
                % output as png graphic via additional matlab function
                if filesave==2
                    savefig([imagefile], 'png', '-rgb', '-c0', '-r250');
                end
            end
        end
    end

disp(['  ']);
disp(['Script ended...']);
end

