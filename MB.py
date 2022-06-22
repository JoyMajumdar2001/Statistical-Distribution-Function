import numpy as np
import matplotlib.pyplot as plt

kb=1.38e-23
e=1.6e-19
E=0.01
dT=0.01

ns=np.arange(3)
n_p=np.arange(100,600,100)
num=int(1000/dT)
T=np.linspace(1,1000,num)
T1=T[:np.size(T)-1].copy()
T2=T1[:np.size(T1)-1].copy()
z=np.empty([np.size(T),np.size(n_p)])
p=np.empty([np.size(T),np.size(ns)])
U=np.empty([np.size(T)-1,np.size(n_p)])
Cv=np.empty([np.size(T)-2,np.size(n_p)])
F=np.empty([np.size(T),np.size(n_p)])
S=np.empty([np.size(T)-1,np.size(n_p)])
for i in range(np.size(n_p)):
    for j in range(np.size(T)):
        zt=0
        for k in range(np.size(ns)):
            p[j,k]=np.exp(-(k*E*e)/(kb*T[j]))
            zt=zt+p[j,k]
        z[j,i]=(zt)**(n_p[i])
        p[j,:]=p[j,:]/zt
    U[:,i]=kb*T1*T1*np.diff(np.log(z[:,i]))/dT
    Cv[:,i]=np.diff(U[:,i])/dT
    F[:,i]=-kb*T*(np.log(z[:,i]))
    S[:,i]=-np.diff(F[:,i])/dT

plt.subplot(2,3,1)  
plt.plot(T,p)
plt.xlabel('Temperature')
plt.ylabel('Probability')
plt.subplot(2,3,2)  
plt.plot(T,z)
plt.xlabel('Temperature')
plt.ylabel('Partition Function')
plt.subplot(2,3,3)  
plt.plot(T1,U)
plt.xlabel('Temperature')
plt.ylabel('Energy')
plt.subplot(2,3,4)  
plt.plot(T2,Cv)
plt.xlabel('Temperature')
plt.ylabel('Specific Heat')
plt.subplot(2,3,5)  
plt.plot(T,F)
plt.xlabel('Temperature')
plt.ylabel('Helmholtz Free Energy')
plt.subplot(2,3,6)  
plt.plot(T1,S)
plt.xlabel('Temperature')
plt.ylabel('Entropy')
plt.show()
