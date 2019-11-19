#!/usr/bin/env python
# coding: utf-8

# In[94]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import subprocess
import os.path


# In[114]:


def generate_data_sho(fName,L,eps,omega,m,hbar):
    result = subprocess.run(['./'+fName,str(L),str(eps),str(omega),str(m),str(hbar)], stdout=subprocess.PIPE)
    splits = result.stdout.decode('utf-8').lstrip().split('\n')
    return (splits[0].strip(), splits[1].strip())


# In[297]:


def do_the_job(L, spacings, omega, m, hbar):
    N=2*L+1
    energies=[]
    states=[]
    x0 = np.sqrt(hbar/(m*omega))
    for spc in spacings:
        enName, stName = generate_data_sho('a.out', L, spc, omega, m, hbar)
        energies.append(np.loadtxt(enName, delimiter=',', comments='#'))
        states.append(np.loadtxt(stName, delimiter=',', comments='#'))
        uvCutoff =(hbar*hbar)/(2*m*spc*spc) 
        print('Done spacing:', spc)
        print('Ratio (UV)', x0/spc )
        print('Ratio (IR)', x0/(spc*L) )
        print('-------------------------------')
    en_test = []
    for i, spc in enumerate(spacings):
        en_test.append((energies[i]-hbar*omega*0.5)/(hbar*omega))
    
    fromN = 1
    toN = N
    fig = plt.figure(figsize=(10,8))
    plt.grid()
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    figTitle = 'Energy vs Energy level L='+str(L)+" $\omega$=" + str(omega)
    plt.title(figTitle)
    plt.plot(np.arange(fromN,toN),np.arange(fromN,toN),c='xkcd:yellow',label='Theoretical', lw=5, ls='--', alpha=0.9)
    for i,ent in enumerate(en_test):
        plt.scatter(np.arange(fromN,toN),ent[fromN:toN], s=10, label='Spacing '+str(spacings[i]))
    plt.legend()
    plt.xlabel('Energy level')
    plt.ylabel('Energy')
    plt.savefig('energyloglog.png')
    plt.show()


# In[298]:


#do_the_job(L=500,spacings = [0.01,0.001, 0.0001],omega=1.0E2,m=1.0,hbar=1.0)


# In[309]:


L=1000
N=2*L+1
spc=0.001
omega=1.0E3
m=1.0
hbar=1.0
enName, stName = generate_data_sho('a.out', L, spc, omega, m, hbar)
E=np.loadtxt(enName, delimiter=',', comments='#')
W=np.loadtxt(stName, delimiter=',', comments='#')


# In[315]:


fig=plt.figure(figsize=(10,8))
plt.grid()
plt.plot(np.arange(N),E, label='Diagonalization', lw=2)
plt.plot(np.arange(N), hbar*omega*(np.arange(N)+0.5), label='Theoretical', lw=2)
plt.legend()
plt.xlabel('Energy level')
plt.ylabel('Energy')
figTitle = 'Energy vs Energy level L='+str(L)+" $\omega$=" + str(omega)
plt.title(figTitle)
plt.savefig('energy.png')
plt.show()


# In[391]:


wfCount = 7
lCut =800
rCut=1200
x=(np.arange(lCut,rCut)-L)*spc
for i in range(wfCount):
    plt.axhline(E[i], c='black', ls='--')
    plt.plot(x, E[i]+3000*W[i*N+lCut:i*N+rCut,0], label='Level '+str(i))
plt.plot(x, 0.5*omega*omega*m*(x**2))
plt.ylim(0, E[wfCount])
plt.xlabel('Position')
plt.ylabel('Potential')
plt.title('First eigenvalues and eigenfunctions (not in scale)')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), shadow=True, ncol=wfCount)
plt.savefig('potential.png')
plt.show()


# In[412]:


wfCount=3
fig,axs = plt.subplots(ncols=1, nrows=wfCount, figsize=(20,10))
for i in range(wfCount):
    psi2 = W[i*N+lCut:i*N+rCut,0]**2+W[i*N+lCut:i*N+rCut,1]**2
    if(np.abs(np.sum(psi2)-1.0)>0.001):
        print('Fatal normalization error')
    axs[i].plot(x,psi2 , label='Level '+str(i), lw=3, c='orangered', ls='-.')
    axs[i].set_ylim(0,0.02)
    axs[i].legend()
fig.suptitle('Probability amplitudes $|\psi|^2$ for the first three levels', fontsize=25)
plt.savefig('pramp.png')
plt.show()


# In[ ]:




