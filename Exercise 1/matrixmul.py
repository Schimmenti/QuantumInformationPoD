import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os.path
sizes = [1,5,10,20,50,100,200,500,750,1000]
times = np.zeros((len(sizes),4))
fName = 'matrixmul.x'
for i,sz in enumerate(sizes):
	result = subprocess.run(['./'+fName,str(sz)], stdout=subprocess.PIPE)
	splits = result.stdout.decode('utf-8').lstrip().split('\n')
	times[i,0]=int(sz)
	times[i,1]=float(splits[0])
	times[i,2]=float(splits[1])
	times[i,3]=float(splits[2])
	print('Matrix dimension is:',sz)
	print('Execution time', times[i,1:])
fig,ax = plt.subplots(nrows=1, ncols=1)
plt.xlabel('Matrix dimension')
plt.ylabel('Log of Execution time [s]')
ax.plot(times[:,0], np.log10(times[:,1]), label='Method 1')
ax.plot(times[:,0], np.log10(times[:,2]), label='Method 2')
ax.plot(times[:,0], np.log10(times[:,3]), label='Intrinsic')
plt.legend()
plt.show()
