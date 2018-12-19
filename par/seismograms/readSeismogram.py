import scipy.io as io
import numpy as np


def read_Seismogram(filename):
    matinfo = io.mminfo('seismogram.p.mtx')
    numrow = matinfo[0]
    numcol = matinfo[1]
    
    mopen = open('seismogram.p.mtx','r')
    readdata = mopen.readlines()
    firstLine = readdata.pop(0)
    secondLine = readdata.pop(0)
    
    plaindata = np.genfromtxt (readdata)
    
    data=np.empty([numrow, numcol])
    
    k=0
    for i in range(numrow):
        for j in range(numcol):
            data[i,j] = plaindata[k]
            k+=1

    return data, numrow, numcol
