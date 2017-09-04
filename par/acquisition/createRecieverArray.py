from writeMtx import write_mtx
import numpy as np

OUTPUT_FILENAME='receiver.mtx'

## Requiered Parameters

## Grid information
NX = 100  # Total number of grid points in X
NY = 100  # Total number of grid points in Y (depth)
NZ = 100  # Total number of grid points in Z

X1 = 20   # X1: First receiver pos in X
X2 = 80   # X2: Last receiber pos in X
Y1 = 20   # Y1: First receiver pos in Y
Y2 = 80   # Y2: Last receiber pos in Y
Z1 = 20   # Z1: First receiver pos in Z
Z2 = 80   # Z2: Last receiber pos in Z
NGEOPH = 1 # Place a geophone on each GP (NGEOPH=1) or every other (NGEOPH=2)

#Define receiver type
RECEIVER_TYPE=[1]   # (1=P, 2=vX, 3=vY, 4=vZ), if you want all write [1 2 3 4]

#Create reciever array
X = np.arange(X1, X2+NGEOPH, NGEOPH)
Y = np.arange(X1, X2+NGEOPH, NGEOPH)
Z = np.arange(X1, X2+NGEOPH, NGEOPH)

len_X = len(X)
RECEIVER_TYPE_vec = np.full((1,len_X),1)[0]

## Write to file

## Create Matrix
RECEIVER_FILE= [X, Y, Z, RECEIVER_TYPE_vec]
print np.shape(RECEIVER_FILE)[0]

## Write mtx file
writeMatrix2mtx(OUTPUT_FILENAME,RECEIVER_FILE);
