function D = merge_snapshots(plot_opts, opts)
%MERGE_SNAPSHOTS  Merge snapshots such that the output has
%   the order of axes (Z, X, Y).
par_folder = plot_opts.par_folder;
file_ext = plot_opts.file_ext;

snap_name = fullfile(par_folder, [opts.SNAP_FILE file_ext]);

nx = opts.NX;
ny = opts.NY;
nz = opts.NZ;

NPROCX = opts.NPROCX;
NPROCY = opts.NPROCY;
NPROCZ = opts.NPROCZ;

IDX = opts.IDX;
IDY = opts.IDY;
IDZ = opts.IDZ;


nlx = (nx/NPROCX)/IDX;
nly = (ny/NPROCY)/IDY;
nlz = (nz/NPROCZ)/IDZ;
nsnap = 1+floor((opts.TSNAP2-opts.TSNAP1)/opts.TSNAPINC);

%%
for i = 1:NPROCX
    disp(i*100/NPROCX);
    for j = 1:NPROCY
        for k = 1:NPROCZ
            fid=fopen([snap_name,'.',num2str(i-1),'.',num2str(j-1),'.',num2str(k-1)]);
            A = fread(fid,'float');
            A = reshape(A,[nly,nlx,nlz,nsnap]);
            A = permute(A,[2,1,3,4]);
            if exist('B','var')
                B = cat(3,B,A);
            else
                B = A;
            end
            clear A
        end
        if exist('C','var')
            C = cat(2,C,B);
        else
            C = B;
        end
        clear B
    end
    if exist('D','var')
        D = cat(1,D,C);
    else
        D = C;
    end
    clear C
end
    
% finally we swap to have Z X Y axis order
D = permute(D,[3, 1, 2, 4]);

end

