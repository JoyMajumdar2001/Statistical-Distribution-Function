import numpy as np
import matplotlib.pyplot as plt

kb=1.38e-23
e=1.6e-19
E=np.array([0.0,0.01,0.02])
dT=0.01
n_p=np.array([1,2,3])
num=int(1000/dT)
T=np.linspace(1,1000,num)
T1=T[:np.size(T)-1].copy()
T2=T1[:np.size(T1)-1].copy()
z=np.empty([np.size(T),np.size(n_p)])
U=np.empty([np.size(T)-1,np.size(n_p)])
Cv=np.empty([np.size(T)-2,np.size(n_p)])
F=np.empty([np.size(T),np.size(n_p)])
S=np.empty([np.size(T)-1,np.size(n_p)])

for j in range(np.size(T)):
    b = -e/(kb*T[j])
    z[j,0] = np.exp(E[0]*b) + np.exp(E[1]*b) + np.exp(E[2]*b)
    z[j,1] = np.exp(2*E[0]*b) + np.exp(2*E[1]*b) + np.exp(2*E[2]*b) + np.exp((E[0] +E[1])*b) + np.exp((E[0] +E[2])*b) + np.exp((E[1] +E[2])*b)
    z[j,2] = np.exp(3*E[0]*b) + np.exp(3*E[1]*b) + np.exp(3*E[2]*b) + np.exp((2*E[0] +E[1])*b) + np.exp((E[0] +2*E[1])*b) + np.exp((2*E[1] +E[2])*b) + np.exp((E[1] +2*E[2])*b) + np.exp((2*E[0] +E[2])*b) + np.exp((E[0] + 2*E[2])*b) + np.exp((E[0] + E[1] + E[2])*b)
for i in range(np.size(n_p)):
    U[:,i]=kb*T1*T1*np.diff(np.log(z[:,i]))/dT
    Cv[:,i]=np.diff(U[:,i])/dT
    F[:,i]=-kb*T*(np.log(z[:,i]))
    S[:,i]=-np.diff(F[:,i])/dT


plt.subplot(2,3,1)  
plt.plot(T,z)
plt.xlabel('Temperature')
plt.ylabel('Partition Function')
plt.subplot(2,3,2)  
plt.plot(T1,U)
plt.xlabel('Temperature')
plt.ylabel('Energy')
plt.subplot(2,3,3)  
plt.plot(T2,Cv)
plt.xlabel('Temperature')
plt.ylabel('Specific Heat')
plt.subplot(2,3,4)  
plt.plot(T,F)
plt.xlabel('Temperature')
plt.ylabel('Helmholtz Free Energy')
plt.subplot(2,3,5)  
plt.plot(T1,S)
plt.xlabel('Temperature')
plt.ylabel('Entropy')
plt.show()
