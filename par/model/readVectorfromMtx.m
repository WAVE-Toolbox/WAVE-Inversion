function [vector]=readVectorfromMtx(filename)

fileID = fopen(filename,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE);
vector=fscanf(fileID,'%e',[1 size(1)]);
fclose(fileID);
end