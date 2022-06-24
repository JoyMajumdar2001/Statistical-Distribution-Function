import numpy as np
import matplotlib.pyplot as plt
import math

e = 1.6e-19
u = 0
kb = 1.38e-23
E = np.arange(-0.5,0.5,0.001)
T = np.arange(100,1100,200)
a = -1
q = 0

FN = np.empty([np.size(E),np.size(T)])
C = np.empty([np.size(FN),3])

for n in range(3):
    for j in range(np.size(T)):
        for i in range(np.size(E)):
            FN[i,j] = 1/(math.exp(((E[i]-u)*e)/(kb*T[j]))+a)
    a = a + 1
    plt.subplot(2,2,n+1)  
    plt.plot(E,FN)
    plt.xlabel('Energy')
    plt.ylabel('f(E)')
