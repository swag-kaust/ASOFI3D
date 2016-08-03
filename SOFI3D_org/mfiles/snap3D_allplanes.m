%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---script for the visualization of snapshots gained from the SOFI3D simulation
%---most parameters are as specified in SOFI3D parameter-file, e.g. sofi3D.json
%---Please note : y dentotes the vertical axis!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---start of input parameter definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---model/snapshot dimensions (gridsize and grid spacing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=100; ny=100; nz=100; % basic grid size; ny=vertical
outx=1; outy=1; outz=1; % snap increment in x/y/z direction, outy=vertical
% spatial discretization, it is assumed that dx=dy=dz=dh
dh=54.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---input, output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose number of input files (1 or 2), only for 2D plots
% For the 3D plot (image_switch=2) only file_inp1 is used by default
num_switch=2;

% Input file1 (snapshot file1)
file_inp1='../par/snap/test.bin.div';
% Input file2 (snapshot file2)
file_inp2='../par/snap/test.bin.curl';
% Model file (for single display or contour plot ontop of snapshot)
file_mod='../par/model/test.SOFI3D.rho';

% Output file
% switch for saving snapshots to picture file 1=yes (jpg) 2= yes (png) other=no
filesave=0;
% base name of picture file output, will be expanded by extension jpg/png
file_out='../par/snap/pic_test';

% title strings for each sub-figure
title_inp1='P-wave field (div)';
title_inp2='S-wave field (curl)';
title_mod='Density model';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---varity of switches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch for contour of model overlaying the model or snapshot file
% 1=yes other=no
cont_switch=1;
% number of contours for the contour-plot
numbOFcont=8;

% Choose model or snapshot plot (model-->1; snapshot-->2)
type_switch=2;
% Choose 2D slice or 3D surf plot
image_switch=2; % 1 = 2D; 2 = 3D;
% Choose slice geometry and postion (for 2-D plots)
slice_switch=2; % horizontal(zx)=1; vertical (yx)=2; vertical (yz)=3;
% slice definition, where to slice trough the 3-D space
nx=nx/outx;ny=ny/outy;nz=nz/outz;
zslice=nz/2; % for xy plane
yslice=ny/2; % for xz plane
xslice=nx/2; % for yz plane in grid points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Snapshot definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time increment for snapshots:
TSNAP1=6.6e-3;
TSNAPINC=0.2;
% firts and last snapshot that is considered for displayin
firstframe=1;
lastframe=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---3D definitions: defines two rotating planes (xz, yz plane)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi1=0; % a horizontal plane (x-z plane) is rotated by phi1 with respect to the rotpoint
phi2=90; % a horizontal plane (x-z plane) is rotated by phi2 with respect to the rotpoint
% rotaxis and rotpoint refers to the rotation of the 2D-slices within the 3D volume
% direction rotation axes [0,1,0] rotation of plane around vertical axis
% [1,0,0] rotation of plane around x-axis
rotaxis=[1,0,0];
% defines point for rotation of the 2D slices [x y z] in meter
%  values are defined as difference to the center of the model/snaphot
% in case of rotpoint=[0,0,0], this point is excatly the center of the model/snaphot
rotpoint=[0 0 0];
% defines angles of perspective for the 3-D view
% i.e. it rotates the plot to favorable perspective
viewpoint=[28,18];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---axis limits for 2D and 3D visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% colorbar boundaries for cropping the snapshot value range
% only used if type_switch=2
auto_scaling=2; % 1= automatic determination of boundaries, 2= constant values caxis_value , 3= no scaling
caxis_value_1=1e-12;
caxis_value_2=1e-12; % only used if num_switch=2

% use custom axis limits if axisoverwrite==1, otherwise matlab will
% determine axis limits automatically
axisoverwrite=0;
newaxis=[35 65 35 65];

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
    num_switch=1;
end

% if 3D snaphot visualization is chosen, the number of snapshots for
% simultaneous display is set to 1
if image_switch==2
    num_switch=1;
end

% load model file ; Format ieee-be for Sun or Workstation, -le for PC
if (cont_switch==1) | (type_switch==1)
    
    % opening file and reading
    disp(['Loading file ' file_mod]);
    fid_mod=fopen(file_mod,'r','ieee-le');
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

if num_switch==2;
    % open snapshot data of 2nd input file
    disp(['Loading snap shot file ' file_inp2]);
    fid_file2=fopen(file_inp2,'r','ieee-le');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---2D display of snapshot data (Slices )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['   ']);

if image_switch==1
    % in case of snapshot files use seismic colormap
    if type_switch==2
        colormap(load('./seismic.map'));
    end
    % creating variables for snapshot content
    % file1_data (and file2_data) depending on number of snapshots
    % displayed simultaneously (switch num_switch)
    
    % selected slice of y-x-plane (vertical plane)
    if(slice_switch==2)
        % allocate memory for 1st snapshot file
        file1_data=zeros(ny,nx);
        if num_switch==2
            % allocate memory for 2nd snapshot file
            file2_data=zeros(ny,nx);
        end
    end
    % selected slice of z-x-plane (horizontal plane)
    if(slice_switch==1)
        % allocate memory for 1st snapshot file
        file1_data2=zeros(nz,nx);
        if num_switch==2
            % allocate memory for 2nd snapshot file
            file2_data2=zeros(nz,nx);
        end
    end
    % selected slice of y-z-plane (vertical plane)
    if(slice_switch==3)
        % allocate memory for 1st snapshot file
        file1_data2=zeros(ny,nz);
        if num_switch==2
            % allocate memory for 2nd snapshot file
            file2_data2=zeros(ny,nz);
        end
    end
    
    % determing imaging vectors (for contour and imagesc) and labels
    % according to chosen image plane
    if slice_switch==1
        image_vect1=x;
        image_vect2=z;
        labelstring1='x in m';
        labelstring2='z in m (horizontal)';
    end
    if slice_switch==2
        image_vect1=x;
        image_vect2=y;
        labelstring1='x in m';
        labelstring2='y in m (vertical)';
    end
    if slice_switch==3
        image_vect1=z;
        image_vect2=y;
        labelstring1='z in m';
        labelstring2='y in m (vertical)';
    end
    
    
    h1=figure(1);
    % determination of screen size
    scrsz = get(0,'ScreenSize');
    % determination size of window for plotting
    set(h1,'Position',[1 scrsz(4)*2/3 scrsz(3)*1/4 scrsz(4)*2/3]);
    
    %creating subfigure handles if 2 snapshots are displayed simultaneously
    if num_switch==2
        ax1=subplot('Position',[0.05 0.3 0.4 0.4]);
        ax2=subplot('Position',[0.55 0.3 0.4 0.4]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---loop over timesteps (2D snapshots)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=firstframe:1:lastframe;
        
        %calculating time of snapshot
        tsnap=(i-1)*TSNAPINC+TSNAP1;
        disp(['Loading snapshot no ',int2str(i),' at time=',num2str(tsnap),' s.']);
     
        % loading data:
        
        if(slice_switch==2) % y-x-plane (vertical)
            % since the models are stores as a series of yx-slices
            % we just have to seek/jump to the mid-slice and load the data
            offset=4*nx*ny*(zslice-1)+4*nx*ny*nz*(i-1);
            fseek(fid_file1,offset,-1);
            if num_switch==2
                fseek(fid_file2,offset,-1);
            end
            file1_data(:,:)=fread(fid_file1,[ny,nx],'float');
            if num_switch==2
                file2_data(:,:)=fread(fid_file2,[ny,nx],'float');
            end
        else
            % z-x-plane (horizontal) and %y-z-plane (vertical)
            % to display slices in any other plane besides y-x, we have to load
            % each single yx slice and extract a single line and put them together
            
            for l=1:ny;
                file1_data=fread(fid_file1,[ny,nx],'float');
                if num_switch==2
                    file2_data=fread(fid_file2,[ny,nx],'float');
                end
                
                if(slice_switch==1) % z-x plane (horizontal)
                    file1_data2(l,:)=file1_data(yslice,:);
                    if num_switch==2
                        file2_data2(l,:)=file2_data(yslice,:);
                    end
                end
                
                if(slice_switch==3) % y-z-plane (vertical)
                    file1_data2(:,l)=file1_data(:,xslice);
                    if num_switch==2
                        file2_data2(:,l)=file2_data(:,xslice);
                    end
                end
                
            end
            file1_data=file1_data2;
            clear file1_data2;
            if num_switch==2
                file2_data=file2_data2;
                clear file2_data2;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---plotting 2D slices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % switch for loading 1 or 2 input files
        if num_switch==2
            % now switching to secondary subplot
            axes(ax2);
            imagesc(image_vect1,image_vect2,file2_data);
            % determing maximum amplitude of plot
            file2max=max(max(abs(file2_data)));
            % switch whether model contour should be plotted
            if cont_switch==1
                hold on
                contour(image_vect1,image_vect2,mod_data,numbOFcont,'k-','LineWidth',1);
                hold off
            end
            
            % formatting the 2nd sub-figure
            colorbar
            xlabel(labelstring1);
            ylabel(labelstring2);
            title(title_inp2);
            set(gca,'DataAspectRatio',[1 1 1]);
            set(get(gca,'title'),'FontSize',12,'FontWeight','bold');
            set(get(gca,'Ylabel'),'FontSize',12,'FontWeight','bold');
            set(get(gca,'Xlabel'),'FontSize',12,'FontWeight','bold');
            set(gca,'FontSize',12,'FontWeight','bold');
            set(gca,'Linewidth',1.0);
            set(gca,'Box','on');
            if axisoverwrite==1
                axis(newaxis);
            end
            
            % limiting the colorbar to specified range
            switch auto_scaling
                case 1
                    caxis([-file2max/10 file2max/10]);
                case 2
                    caxis([-caxis_value_2 caxis_value_2])
                otherwise
            end
            
            % now switching to primary subplot
            axes(ax1);
        end
        
        % plot 1st input file
        if type_switch==2
            % plot snapshot file
            imagesc(image_vect1,image_vect2,file1_data);
        end
        if type_switch==1
            % plot model file
            imagesc(image_vect1,image_vect2,mod_data);
        end
        colorbar
        
        if type_switch==2
            title(title_inp1);
        end
        if type_switch==1
            title(title_mod);
        end
        % determing maximum amplitude of plot
        file1max=max(max(abs(file1_data)));
        
        % switch whether model contour should be plotted
        if cont_switch==1
            hold on
            contour(image_vect1,image_vect2,mod_data,numbOFcont,'k-','LineWidth',1);
            hold off
        end
        
        % formating the 1st sub-figure
        xlabel(labelstring1);
        ylabel(labelstring2);
        set(gca,'DataAspectRatio',[1 1 1]);
        set(get(gca,'title'),'FontSize',12,'FontWeight','bold');
        set(get(gca,'Ylabel'),'FontSize',12,'FontWeight','bold');
        set(get(gca,'Xlabel'),'FontSize',12,'FontWeight','bold');
        set(gca,'FontSize',12,'FontWeight','bold');
        set(gca,'Linewidth',1.0);
        set(gca,'Box','on');
        if axisoverwrite==1
            axis(newaxis);
        end
        
        if type_switch==2 % in case of snapshot:
            % display maximum amplitude of sub-figures to command window
            if num_switch==2
                disp(['  Maximum amplitude of ',title_inp2,'-snapshot: ', num2str(file2max)]);
            end
            disp(['  Maximum amplitude of ',title_inp1,'-snapshot: ', num2str(file1max)]);
            
            % limiting the colorbar to specified range
            switch auto_scaling
                case 1
                    caxis([-file1max/10 file1max/10]);
                case 2
                    caxis([-caxis_value_1 caxis_value_1]);
                otherwise
            end
            
            % adding text string for timestep of snapshot
            set(text(7000,7000,['T2= ',sprintf('%2.4f',tsnap), 's']),'FontSize',14,'FontWeight','bold','Color','black');
        end
        
        % delay the display of snapshots by frictions of seconds
        pause(pause4display);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Saving the snapshot to file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (filesave~=0)
            % generating filename string
            if i<10
                imagefile=[file_out,'2D_00',int2str(i)];
            else if i<100
                    imagefile=[file_out,'2D_0',int2str(i)];
                else
                    imagefile=[file_out,'2D_',int2str(i)];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---3D display of snapshot data (single snapshot only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if image_switch==2
    
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
        colormap(load('./seismic.map'));
        
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
    for i=firstframe:1:lastframe;
        
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
        % creating a grid
        [X,Z,Y]=meshgrid(x,z,y);
        
        %surface height of horizontal z-x plane
        yplane=zeros(ny,nx);
        yplane(:,:)=(ny*dh*outy)/2+rotpoint(2); 
        
        %surface height of vertical y-x plane
        zplane=zeros(nz,nx);
        zplane(:,:)=(nz*dh*outz)/2+rotpoint(3); 
        
        % creating a vertical slice plane in y-x plane
        hslice = surf(x,y,yplane);
        % rotate slice plane with, with respect to a point from which rotation is defined
        % in case of rotpoint=[0,0,0], this point is the center of the model/snaphot
        rotate(hslice,rotaxis,phi1,[mean(x)-rotpoint(1),mean(y)-rotpoint(2),mean(mean(zplane))-rotpoint(3)]);
        % get boundaries of rotated plane
        xd = get(hslice,'XData');
        yd = get(hslice,'YData');
        zd = get(hslice,'ZData');
        % remove plane, only limits are further used
        delete(hslice);
        
        % creating a horizontal slice plane in z-x plane
        hslice2 = surf(x,z,zplane);
        % rotate slice plane with, with respect to a point from which rotation is defined
        % in case of rotpoint=[0,0,0], this point is the center of the model/snaphot
        rotate(hslice2,rotaxis,phi2,[mean(x)-rotpoint(1),mean(z)-rotpoint(2),mean(mean(yplane))-rotpoint(3)]);
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
        set(h1,'Position',[1 scrsz(4)*2/3 scrsz(3)*1/4 scrsz(4)*2/3]);
        
        % strict vertical slice, no rotation applied
%                 h = slice(Y,X,Z,file1_data,[],yslice*dh*outy-1,[]); 
        % rotated vertical slice
        % !!! note that taking sclices from a homogeneous model is - for some
        % reason not working, please use 2-D visualization instead !!!
        h = slice(X,Z,Y,file1_data,xd,zd,yd);
        set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
        
        hold on
        % strict horizontal, no rotation applied
        %         h2 = slice(Y,X,Z,file1_data,[],[],zslice*dh*outz); % this is strict vertical, no rotation applied
        % rotated horizontal slice
        h2 = slice(X,Z,Y,file1_data,xd2,zd2,yd2);
        set(h2,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
        
        % vertical slice at model boundary
        %         h3 = slice(X,Z,Y,file1_data,10,[],[]);
        %         set(h3,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
        
        % black outline of the vertical slice
        plot3([0 max(max(xd))],[max(max(zd)) max(max(zd))],[0 0],'-black','LineWidth',2);
        plot3([max(max(xd)) max(max(xd))],[max(max(zd)) max(max(zd))],[0 max(max(yd))],'-black','LineWidth',2);
        plot3([max(max(xd)) 0],[max(max(zd)) max(max(zd))],[max(max(yd)) max(max(yd))],'-black','LineWidth',2);
        plot3([0 0],[max(max(zd)) max(max(zd))],[0 max(max(yd))],'-black','LineWidth',2);
        
        % black outline of the horizontal slice
        plot3([0 max(max(xd2))],[max(max(zd2)) max(max(zd2))],[max(max(yd2)) min(min(yd2))],'-black','LineWidth',2);
        plot3([max(max(xd2)) max(max(xd2))],[max(max(zd2)) 0],[max(max(yd2)) max(max(yd2))],'-black','LineWidth',2);
        plot3([max(max(xd2)) 0],[0 0],[max(max(yd2)) min(min(yd2))],'-black','LineWidth',2);
        plot3([0 0],[0 max(max(zd2))],[max(max(yd2)) max(max(yd2))],'-black','LineWidth',2);
        
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
        daspect([1,1,1]);
        axis tight
        box on
        
        % set viewpoint within 3-D space
        % (rotate 3-D figure plot to favorable perspective)
        view(viewpoint);
        camproj perspective
        % lightangle(-45,45);
        
        set(gcf,'Renderer','zbuffer');
        colorbar
        xlabel('x in m');
        ylabel('z in m');
        zlabel('Depth y in m');
        
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
                    caxis([-file1max/10 file1max/10]);
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
end
disp(['  ']);
disp(['Script ended...']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---End of Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%