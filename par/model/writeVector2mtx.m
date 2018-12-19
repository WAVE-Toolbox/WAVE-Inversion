function writeVector2mtx(filename,vector)

n=size(vector,1);

fileID = fopen(filename,'w');
fprintf(fileID,['%%']);
fprintf(fileID,['%%MatrixMarket vector array real general\n']);
fprintf(fileID,[num2str(numel(vector)), ' ', num2str(1)]);
fprintf(fileID,'\n');
fprintf(fileID,'%f \n',[vector]);
fclose(fileID);

end