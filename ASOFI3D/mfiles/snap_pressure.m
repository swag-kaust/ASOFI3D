close all
clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
filediv='test_Khang_cut_plane.bin.p';

% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=100;
NY1=1; NY2=100; 
IDX=1; IDY=1;
dh=54.0;

% time increment for snapshots:
TSNAPINC=0.2; TSNAP1=6.6e-3;
FW=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

firstframe=1;
lastframe=25;
snapno=0;

caxis_value=1.0e-3;

load 'seismic.map'
colormap(seismic)

% load model
% file='adapt_test_Khang_model.bin';
% disp([' loading file ' file]);
% fid=fopen(file,'r','ieee-le');
% model=fread(fid,[ny,nx],'float');
% fclose(fid);
 %model=model(1:2:size(model,1),1:2:size(model,2));
 
 disp(['opening file ' filediv]);
 fid_div=fopen(filediv,'r','ieee-le');

 for i=firstframe:1:lastframe,
 
 hold off
   disp(['loading snapshot no ',int2str(i)]);
   % loading data:
    tsnap=(i-1)*TSNAPINC+TSNAP1;
   offset=4*nx*ny*(i-1);
   fseek(fid_div,offset,-1);
   veldiv=fread(fid_div,[ny,nx],'float');
   vmp=max(max(abs(veldiv)));
   disp([' Maximum amplitude of P-snapshots: ', num2str(vmp)]);

		 
	set(gcf,'Color',[1 1 1]) 
       imagesc(x,y,veldiv);  
		 hold on
 		hold off
		 if (i>0), caxis([-caxis_value caxis_value]); end
		
       set(text(-1500.0,0.0,['T=',sprintf('%1.2f',tsnap*1000.0),' ms']),...
	  'FontSize',12,'FontWeight','bold','color','k');  
       title('P-Waves')
       xlabel('Distance [m]')
       ylabel('Depth [m]')
       set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
 

   	pause(0.1);
	

	snapno=snapno+1;
	
    % Saving the snapshot:
    pngfile=['Kang_',int2str(i),'.png'];
    eval(['print -dpng ' pngfile]);
    
 end
fclose(fid_div);











