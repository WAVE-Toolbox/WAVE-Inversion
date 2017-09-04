from readSeismogram import read_Seismogram
import matplotlib.pyplot as plt
import numpy as np

#Read Seismogram
filename = 'seismogram.p.mtx'
data, numrow, numcol =read_Seismogram(filename)

#Plot Seismogram
plt.figure(facecolor='white')
x = np.arange(0,numrow,1)
filling = False

for i in range(numcol):
    yi = data[:,i]
    normalized_yi = yi / abs(yi.max(axis=0)) +i+1
    if filling:
        plt.fill_between(x, i+1, normalized_yi, facecolor='black')
        plt.plot(x,normalized_yi, c='k', linewidth=.5)
    else:
        plt.plot(x,normalized_yi, c='k')

#Plot-Style
plt.title('Normalized traces')
plt.xlabel('Samples')
plt.ylabel('Traces')
plt.xlim([0,numrow+1])
plt.ylim([0,numcol])

plt.show()

