#!/usr/bin/env python
# coding: utf-8

# In[151]:


import numpy as np
import matplotlib.pyplot as plt


# In[167]:


Nmax = 10
spectra = []
for i in range(2,Nmax+1):
    spectra.append(np.loadtxt('isingsp_'+str(i)+'.txt', delimiter=' '))
    #spectra.append(np.loadtxt('anderson_250.txt', delimiter=' '))


# In[170]:


fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.linspace(0,2), -1-0.25*np.linspace(0,2)**2, c='black')
ax.plot(np.linspace(2,3),-np.linspace(2,3), c='black', label='Mean field')
ax.set_xlabel('$\lambda$')
ax.set_ylabel('Ground state/N-1')
for i, sp in enumerate(spectra):
    lmbd = sp[:,0]
    eigs = sp[:,1:]
    plt.plot(lmbd, eigs[:,0]/(i+1), label='N='+str(i+2))
plt.legend()
plt.grid()
plt.title('Ground states')
plt.savefig('images/groundstates.png')
plt.show()


# In[233]:


for i, sp in enumerate(spectra):
    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_xlabel('$\lambda$')
    ax.set_ylabel('Energy/N-1')
    lmbd = sp[:,0]
    eigs = sp[:,1:]
    for j in range(min(eigs.shape[1],4)):
        plt.plot(lmbd, eigs[:,j]/(i+1), label='$e_'+str(j)+"$")
    plt.title('First four eigenvalues for $N='+str(i+2)+'$')
    plt.legend()
    plt.grid()
    if(i != 0):
        plt.ylim(-3.5, 0.5)
    plt.savefig('images/eigenvalues'+str(i+2)+'.png')
    plt.show()


# In[158]:


def buildHamiltonian(N):
    H = np.zeros((2**N, 2**N))
    for i in range(2**N):
        mask = 3
        for k in range(N-1):
            H[i,i^mask]=1
            H[i^mask,i]=1
            mask = mask*2
    return H


# In[186]:


for nn in range(2,14):
    fig, ax = plt.subplots(figsize=(10,8))
    ax.matshow(buildHamiltonian(nn))
    plt.title('Hamiltonian matrix elements for $N='+str(nn)+'$', y=1.08)
    plt.savefig('images/mat'+str(nn)+'.png')
    plt.show()


# In[205]:


Nmax = 10
gstates = []
for i in range(2,Nmax+1):
    gstates.append(np.loadtxt('isinggs_'+str(i)+'.txt', delimiter=' '))

