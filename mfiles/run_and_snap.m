function run_and_snap


%% run modeling
cd ../
system('./run_asofi3D.sh 96')
cd mfiles

%% show snapshots

[opts, plot_opts] = snap3D_ASOFI;


end
