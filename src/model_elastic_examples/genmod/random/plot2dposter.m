% parameter
NX=240; NY=600; DH=2.0;
file='random2d2.vp';
%--------------------------

close all
fp=fopen(file,'r');
vp=fread(fp,[NY NX],'float');
imagesc(vp);
set(gca,'DataAspectRatio',[1 1 1]); 	
axis off
colormap(gray)
% print -depsc -r300 -zbuffer random2d.eps 	
	
