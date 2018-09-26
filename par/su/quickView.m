close all
clear all

% the script reads the seismogram (single shot) and takes a slice of its 3D Fourier 
% the central frequency is shifted towards the upper limit, because of the
% ghost ?!?
%% parameters

% spatial and time sampling
dx = 50; % meters
dt = 5*10^-3; % seconds
mp = 17; % source is located at mp,mp grid point
mute_p = 30; % mute in the square around mp
azimuth = 0;

%% 

fid = fopen('test_p.bin');
a = fread(fid,'float');
sz_hor = sqrt(length(a)/200);
a = reshape(a,[200 sz_hor sz_hor]);
%a = a(:,mp:end,mp:end);
%a(:,mp-mute_p:mp+mute_p,mp-mute_p:mp+mute_p) = 0;
fclose(fid);
n = size(a);


% prepare axis for fourier domain
fSpace.t = linspace(0,dt*n(1),n(1));
fSpace.x = linspace(-dx*mp,dx*(n(2)-mp),n(2));

% geometrical spreading correction * |distance|

figure
%close all;
subplot(1,2,1)
imagesc(fSpace.x,fSpace.t,squeeze(a(:,mp,:)));
cax = caxis();
caxis(cax/10000);
title('original data')
colormap(load('./srgb.map'));


