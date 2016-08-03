% parameter
NX=240; NY=600; DH=2.0;
file='random2d2.vp';
titel='a=45m, nu=0.5, sigma=5 Prozent';
%--------------------------

close all

x=(1:NX)*DH;
y=(1:NY)*DH;

fp=fopen(file,'r');
vp=fread(fp,[NY NX],'float');
imagesc(x,y,vp);
xlabel(' X[m]');
ylabel(' Y[m]');
title(titel);
set(gca,'DataAspectRatio',[1 1 1]); 	set(get(gca,'title'),'FontSize',12);
 	
	
	set(get(gca,'title'),'FontWeight','bold');
 	set(get(gca,'Ylabel'),'FontSize',16);
 	set(get(gca,'Ylabel'),'FontWeight','bold');
 	set(get(gca,'Xlabel'),'FontSize',16);
 	set(get(gca,'Xlabel'),'FontWeight','bold');
 	set(gca,'Box','on');
 	set(gca,'FontSize',16);
 	set(gca,'FontWeight','bold');
 	set(gca,'Linewidth',1.0);
