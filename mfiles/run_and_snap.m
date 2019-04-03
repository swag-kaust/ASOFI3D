function run_and_snap


%% run modeling
% load config
config = read_asofi3D_json('../par/in_and_out/asofi3D_hom.json');
jsonPath = '../par/in_and_out/asofi3D.json';

config.jsonPath = jsonPath;
config_hom = config;

%%  generate homogeneous data
config.DH1 = 1e6;
write_asofi3D_json(jsonPath, config);

% run modeling
run_asofi3D(config_hom);



%% show snapshots
[opts, plot_opts, D] = snap3D_ASOFI('hom');

save('D_hom','D','-v7.3');

%%

% isotropic parameters are perturbed by 5%
iso_par = {'VPV2', 'VSV2', 'RHO2'};
for i = 1:length(iso_par)
    config = config_hom;
    config.(iso_par{i}) = 1.1 * config.(iso_par{i});
    config.folderOut = iso_par{i};
    last_snap = create_snapshots(config);
    last_snap_collection(i, :) = last_snap(:);
end

%% anisotropic deviation parameters are perturbed by 0.1
aniso_par = {'EPSX2', 'EPSY2', 'DELX2', 'DELY2', 'DELXY2', 'GAMX2', 'GAMY2'};
for i = 1:length(aniso_par)
    config = config_hom;
    % perturb parameter
    config.(aniso_par{i}) = 0.1;
    % saving figure to folderOut
    config.folderOut = aniso_par{i};
    last_snap = create_snapshots(config);
    last_snap_collection(i+3, :) = last_snap(:);
end

%% SVD

[U, S, V] = svd(last_snap_collection(:,:),'econ');
figure(777);
S = diag(S);
semilogy(S/max(S));
ylim([10^-2 1]);
%semilogy((S(1:end-1)-S(2:end))/max(S));
%imagesc(U);

%% full model
write_asofi3D_json(jsonPath, config_hom);

end

%%
function run_asofi3D(config)
% write config to jsonPath

write_asofi3D_json(config.jsonPath, config);


%% run modeling
% find number of processors NP to run on
NP = config.NPROCX * config.NPROCY * config.NPROCZ;
cd ../
system(['./run_asofi3D.sh ', num2str(NP)])
cd mfiles

end

function last_snap = create_snapshots(config)
% runs asofi3D, then generates snapshots

run_asofi3D(config);

[opts, plot_opts, D] = snap3D_ASOFI(config.folderOut, true);
caxis(caxis());

save('D_diff','D','-v7.3');

%%

last_snap = squeeze(D(:,:,:,3));


[M, I] = max(abs(last_snap), [], 3);


[ind_1, ind_2] = ndgrid(1:size(I,1), 1:size(I,2));



I_linear = sub2ind(size(last_snap), ind_1(:), ind_2(:), I(:));
figure; surf(I, reshape(last_snap(I_linear),size(I)));
set(gca,'Zdir','reverse')
colormap(rdbuMap());
title('Max along Z');

%%


M = M(:);
mask_wave = M > max(M)/10;

K = abs(last_snap(I_linear));
figure; scatter3(ind_1(mask_wave), ind_2(mask_wave), I(mask_wave), (K(mask_wave)./max(K))+eps);
set(gca,'Zdir','reverse')
axis equal
xlim([0 size(last_snap,1)])
ylim([0 size(last_snap,2)])
zlim([0 size(last_snap,3)])

%% 3D scatter snapshot visualization
last_snap = squeeze(D(:,:,:,3));




[Y, X, Z] = ndgrid(0:size(last_snap,1)-1, ...
    0:size(last_snap,2)-1, ...
    0:size(last_snap,3)-1);

X = config.IDX * config.DX * X;
Y = config.IDY * config.DY * Y;
Z = config.IDZ * config.DZ * Z;

mask_3D = (abs(last_snap(:)) > max(abs(last_snap(:)))*0.2);% .* (Z(:)<config.DH1);

mask_3D = logical(mask_3D);

figure; scatter3(X(mask_3D), Y(mask_3D), Z(mask_3D), 100*abs(last_snap(mask_3D))./max(last_snap(:))+eps, last_snap(mask_3D), 'filled');
axis equal
xlim([0 max(X(:))])
xlabel X
ylim([0 max(Y(:))])
ylabel Y
zlim([0 max(Z(:))])
zlabel Z
colormap(rdbuMap());
colorbar
title(config.folderOut);

caxis(caxis()/10);
set(gca,'Zdir','reverse')

end








