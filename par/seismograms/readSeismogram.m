
function [Seismogram]=read_seismogram(filename)

fileID = fopen(filename,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE);
A=fscanf(fileID,'%e',[1 size(1)*size(2)]);
Seismogram=reshape(A(:),[size(1), size(2)]);

end