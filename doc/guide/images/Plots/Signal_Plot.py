import matplotlib.pyplot as plt
import numpy as np
from matplotlib import interactive
interactive(True)

t= np.arange(0.,100,0.01)
FC = 0.05
Tshift = 20
Amp = 1

tau1 = np.pi*FC*(t - 1.2/FC)
f1=-Amp * 2 * tau1 * np.exp(-tau1 * tau1)
maxf1=max(f1)

t2_0= np.arange(0.,20,0.01)
t2_1= np.arange(20.,40,0.01)
t2_2= np.arange(40.,100,0.01)
tau2 = np.pi* 2 *FC*(t2_1 - Tshift)
f2_1 = Amp * (0.5 - 0.75 * np.cos(tau2) + 0.25 * (np.cos(tau2))*(np.cos(tau2))*(np.cos(tau2)))/FC/0.75/np.pi
f2_2=0*t2_2
f2_0=0*t2_0
maxf2 = max(f2_1)

tau3 = np.pi*FC*(t - 1.5/FC)
f3 = Amp*(1 - 2*tau3*tau3)*np.exp(-tau3*tau3)

t4_1= np.arange(20.,40,0.01)
t4_2= np.arange(40.,100,0.01)
t4_3= np.arange(0.,20,0.01)
tau4 = np.pi*FC*(t4_1 - Tshift)
f4_1 = Amp* np.sin(tau4)* np.sin(tau4)* np.sin(tau4)
f4_2=0*t4_2
f4_3=0*t4_3

t5_1= np.arange(20.,40,0.01)
t5_2= np.arange(40.,100,0.01)
t5_3= np.arange(0.,20,0.01)
tau5 = 2*np.pi*FC*(t5_1 - Tshift)
f5_1 = Amp*(np.sin(tau5) - 0.5 * np.sin(2*tau5))
maxf5=max(f5_1)
f5_2=0*t5_2
f5_3=0*t5_3

f6=0*t

plt.plot( t, f3, '-g', lw=1, label="Ricker")

plt.plot( t5_1, f5_1/maxf5, '-b', lw=1, label="Fuchs-Mueller")
plt.plot( t5_2, f5_2, '-b', lw=1)
plt.plot( t5_3, f5_3, '-b', lw=1)

plt.plot( t4_1, f4_1, ':k', lw=1, label="SinThree")
plt.plot( t4_2, f4_2, ':k', lw=1)
plt.plot( t4_3, f4_3, ':k', lw=1)

plt.plot(t, f1/maxf1, '--r', dashes=(2, 1), lw=1,  label="FGaussian")

plt.vlines(x=20, ymin=0., ymax =1, lw=1, color = 'm')
plt.plot( t, f6, 'm', lw=1, label="Spike")

plt.plot(t2_1, f2_1/maxf2, dashes=[2, 1, 1, 1, 1, 1] ,color = 'c', lw=1, label="IntSinThree")
plt.plot(t2_2, f2_2, dashes=[2, 1, 1, 1, 1, 1] ,color = 'c', lw=1)
plt.plot(t2_0, f2_0, dashes=[2, 1, 1, 1, 1, 1] ,color = 'c', lw=1)

plt.legend(loc='upper right')
plt.ylabel('Amplitude')
plt.xlabel('Time [ms]')
plt.xlim((0,100))

ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.show()
raw_input('press return to end')
plt.close('all')
