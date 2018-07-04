% Input file1 (snapshot file1)
file_inp1='../par/snap/test.bin.div';
% Input file2 (snapshot file2)
file_inp2='../par/snap/test.bin.curl';
% Model file (for single display or contour plot ontop of snapshot)
file_mod='../par/snap/test.bin.div';

%dimensions of binary
nXYZT = [512 512 512 5];

A = zeros(nXYZT);

fid_mod=fopen(file_inp1,'r','ieee-le');
A=fread(fid_mod,'float');
A=reshape(A,nXYZT);
fclose(fid_mod);

%slice(X,Z,Y,file1_data,xd2,zd2,yd2);

%%
% % tic;
% % nFrac=10;
% % B3 = zeros(nXYZT(1),nXYZT(2),nXYZT(3));
% % 
% % for i=1:nXYZT(1)/4
% %     B3_i = zeros(nXYZT(1),nXYZT(2));
% %     parfor j=1:nXYZT(2)/4
% %         for k=1:nXYZT(3)/4
% %             xq = 0:X(i,j,k)/nFrac:X(i,j,k);
% %             yq = 0:Y(i,j,k)/nFrac:Y(i,j,k);
% %             zq = 0:Z(i,j,k)/nFrac:Z(i,j,k);
% %             B3_i(j,k) = sum(B_fun(xq,yq,zq));
% %         end
% %     end
% %     B3(i,:,:) = B3_i(:,:);
% % end
% % toc;