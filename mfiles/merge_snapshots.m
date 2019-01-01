function D = merge_snapshots(plot_opts, opts)
%MERGE_SNAPSHOTS  Merge snapshots such that the output has
%   the order of axes (Z, X, Y).
par_folder = plot_opts.par_folder;
file_ext = plot_opts.file_ext;

snap_name = fullfile(par_folder, [opts.SNAP_FILE file_ext]);

nlx = (opts.NX / opts.NPROCX) / opts.IDX;
nly = (opts.NY / opts.NPROCY) / opts.IDY;
nlz = (opts.NZ / opts.NPROCZ) / opts.IDZ;

TSNAP1 = opts.TSNAP1;
TSNAP2 = opts.TSNAP2;
TIME = opts.TIME;
TSNAPINC = opts.TSNAPINC;

if TSNAP2 > TIME
    fprintf(['WARNING: TSNAP2 = %f is larger than TIME = %f. ' ...
             'Set TSNAP2 = TIME.\n'], ...
            TSNAP2, TIME);
    TSNAP2 = TIME;
end

nsnap = 1+floor((TSNAP2 - TSNAP1) / TSNAPINC);

%%
for i = 1:opts.NPROCX
    disp(i*100/opts.NPROCX);
    for j = 1:opts.NPROCY
        for k = 1:opts.NPROCZ
            snap_file = [snap_name,'.',num2str(i-1),'.',num2str(j-1),'.',num2str(k-1)];
            fid=fopen(snap_file);
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
