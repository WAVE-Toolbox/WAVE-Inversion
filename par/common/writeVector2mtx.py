def writeVector2mtx(filename, array):

    fileID = open(filename, 'w')
    fileID.write('%%MatrixMarket vector array real general\n')
    fileID.write(str(len(array))+' 1\n')
    fileID.write('\n'.join(map(str, array)))
    fileID.write('\n')
    fileID.close()

