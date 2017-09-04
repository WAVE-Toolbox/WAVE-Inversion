function writeMatrix2mtx(filename,matrix)

n=size(matrix,1);

[row,col,v] = find(matrix);
fileID = fopen(filename,'w');
fprintf(fileID,['%%']);
fprintf(fileID,['%%MatrixMarket matrix coordinate real general\n']);
fprintf(fileID,[num2str(size(matrix,1)), ' ', num2str(size(matrix,2)), ' ', num2str(numel(v))]);
fprintf(fileID,'\n');
if(n>1)
fprintf(fileID,'%i %i %i\n',[row, col, v]');
else
fprintf(fileID,'%i %i %i\n',[row; col; v]);
end
fclose(fileID);

end