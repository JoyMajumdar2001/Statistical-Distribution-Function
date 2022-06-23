import numpy as np
import matplotlib.pyplot as plt

h = 6.626e-34
kb = 1.38e-23
c = 3e8
T = np.arange(500,1000,5)
T0 = 1000
Wn = np.arange(100,30000,5)
W = Wn*10**(-9)
A = 8*3.14*h*c
q=0

UP = np.empty([np.size(W),np.size(T)])
UR = np.empty([np.size(W),np.size(T)])
UW = np.empty([np.size(W),np.size(T)])

for j in range(np.size(T)):
    for i in range(np.size(W)):
        UP[i,j] = A/(W[i]**(5))/(np.exp(h*c/(W[i]*kb*T[j])) - 1)
        UR[i,j] = 8*3.14*kb*T[j]/(W[i]**4)
        UW[i,j] = (A/(W[i]**(5)))*(np.exp(-1*h*c/(W[i]*kb*T[j])))
    if T[j]==T0:
        q=j
        
plt.subplot(2,2,1)  
plt.plot(Wn,UP)
plt.xlabel('Wavelength')
plt.ylabel('Energy(Plank)')
plt.subplot(2,2,2)  
plt.plot(Wn,UR)
plt.ylim([0, 15])
plt.xlabel('Wavelength')
plt.ylabel('Energy(Rayleigh)')
plt.subplot(2,2,3)  
plt.plot(Wn,UW)
plt.xlabel('Wavelength')
plt.ylabel('Energy(Weins)')
plt.subplot(2,2,4)  
plt.plot(Wn,UP[:,q])
plt.plot(Wn,UR[:,q])
plt.plot(Wn,UW[:,q])
plt.ylim([0, 7])
plt.xlabel('Wavelength')
plt.ylabel('Energy')
plt.show()
