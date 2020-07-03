import numpy as np
import matplotlib.pyplot as plt


X = np.loadtxt("gstate_4_1000.txt",delimiter=',')
fig,ax=plt.subplots(figsize=(10,8))
ax.plot(X[:,0], X[:,1])
ax.set_xlabel("$\lambda$")
ax.set_ylabel('Ground state density')
ax.grid()
plt.savefig('gstate.png')