import numpy as np

def writeMatrix2mtx(filename, matrix):

    n = np.shape(matrix)[1]
    m = np.shape(matrix)[0]

    fileID = open(filename, 'w')
    fileID.write('%%MatrixMarket vector array real general\n')
    fileID.write(str(n) + ' ' + str(m) + ' ' + str(np.count_nonzero(matrix)) + '\n')
    if n > 1:
        for i in range(m):
            for j in range(n):
                if matrix[i][j] != 0:
                    fileID.write(str(j+1) + ' ' + str(i+1) + ' ' + str(matrix[i][j]) + '\n')
    else:
        for i in range(n):
         for j in range(m):
             if matrix[j][0] != 0:
                 fileID.write(str(i+1) + ' ' + str(j+1) + ' ' + str(matrix[j][0]) + '\n')
    FileID.write('\n')
    fileID.close()



#
#function write_mtx(filename,matrix)
#
#n=size(matrix,1);
#
#[row,col,v] = find(matrix);
#fileID = fopen(filename,'w');
#fprintf(fileID,['%%']);
#fprintf(fileID,['%%MatrixMarket matrix coordinate real general\n']);
#fprintf(fileID,[num2str(size(matrix,1)), ' ', num2str(size(matrix,2)), ' ', num2str(numel(v))]);
#fprintf(fileID,'\n');
#if(n>1)
#fprintf(fileID,'%i %i %i\n',[row, col, v]');
#        else
#        fprintf(fileID,'%i %i %i\n',[row; col; v]);
#        end
#        fclose(fileID);
#        
#        end
