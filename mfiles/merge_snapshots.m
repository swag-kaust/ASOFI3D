%% merges snapshots otput has Z X Y order of axis
function D = merge_snapshots(par_folder)

%%
% clear all
% par_folder = '../par/';
par_folder = [par_folder, '/'];
jV = read_asofi3D_json([par_folder, 'in_and_out/sofi3D.json']);

snap_name = [par_folder, jV.SNAP_FILE,'.bin.div'];

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

