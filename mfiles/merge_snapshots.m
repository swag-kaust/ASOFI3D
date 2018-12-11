function D = merge_snapshots(plot_opts, opts)
%MERGE_SNAPSHOTS  Merge snapshots such that the output has
%   the order of axes (Z, X, Y).
par_folder = plot_opts.par_folder;
file_ext = plot_opts.file_ext;

snap_name = fullfile(par_folder, [opts.SNAP_FILE file_ext]);

nx = str2num(jV.NX);
ny = str2num(jV.NY);
nz = str2num(jV.NZ);

NPROCX = str2num(jV.NPROCX);
NPROCY = str2num(jV.NPROCY);
NPROCZ = str2num(jV.NPROCZ);

IDX = str2num(jV.IDX);
IDY = str2num(jV.IDY);
IDZ = str2num(jV.IDZ);


nlx = (nx/NPROCX)/IDX;
nly = (ny/NPROCY)/IDY;
nlz = (nz/NPROCZ)/IDZ;
nsnap = 1+floor((str2num(jV.TSNAP2)-str2num(jV.TSNAP1))/str2num(jV.TSNAPINC));

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

