% create colormap

a = zeros(2001,3);
a(:,:)=NaN;
a(1,:) = [1 0 0];
a(1001,:) = [0.85 1 0.85];
a(2001,:) = [0 0 1];
a = fillmissing(a,'linear',1);
close all;
imagesc(a);
dlmwrite('srgb.map',a);


