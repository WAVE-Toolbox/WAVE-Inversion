function [model]=readModelfromMtx(filename,NX,NY,NZ)

fileID = fopen(filename,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE);
A=fscanf(fileID,'%e',[1 size(1)]);
model=permute((reshape(A(:),[NX, NY, NZ])),[2 1 3]);

end