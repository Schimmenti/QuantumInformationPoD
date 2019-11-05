import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

if(len(sys.argv) >= 4):
	N_min = int(sys.argv[1])
	N_max = int(sys.argv[2])
	N_count = int(sys.argv[3])
else:
	N_min=50
	N_max=1500
	N_count = 50
sizes = np.logspace(np.log10(N_min), np.log10(N_max+1),num=N_count, base=10.0, dtype=int)
times = np.zeros((len(sizes),7))
fName = 'matrixmul.x'
modes = ['ijk','ikj','jik','jki','kij','kji']
for i,sz in enumerate(sizes):
	print('%.2f' % float(i*100/len(sizes)),'%')
	print('Matrix size:',sz)
	result = subprocess.run(['./'+fName,str(sz)], stdout=subprocess.PIPE)
	splits = result.stdout.decode('utf-8').lstrip().split('\n')
	times[i,:]=np.array([float(tm) for tm in splits[0:7]])
doPlots = False

for i in range(6):
	np.savetxt(fname='mode'+str(modes[i]),X=np.concatenate((sizes.reshape(-1,1),times[:,i].reshape(-1,1)),axis=1), fmt=['%i','%.3E'])
np.savetxt(fname='intrinsic',X=np.concatenate((sizes.reshape(-1,1),times[:,6].reshape(-1,1)),axis=1), fmt=['%i','%.3E'])
if(doPlots):
	for i in range(6):
		plt.plot(sizes, times[:,i], label='Mode ' + modes[i])
		plt.plot(sizes, times[:,-1], label='Instrinsic')
		plt.xlabel('Dimension')
		plt.ylabel('Execution time [s]')
		plt.title('Execution time of matrix multiplication')
		plt.legend()
		plt.savefig('extime.png')
		plt.show()


	for i in range(6):
		plt.plot(sizes, np.log10(times[:,i]), label='Mode ' + modes[i])
		plt.plot(sizes, np.log10(times[:,-1]), label='Instrinsic')
		plt.xlabel('Dimension')
		plt.ylabel('Log execution time [s]')
		plt.title('Log execution time of matrix multiplication')
		plt.legend()
		plt.savefig('logextime.png')
		plt.show()
