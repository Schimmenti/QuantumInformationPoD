import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

def Ps(x, a,b,c,d):
    return a*(x**b)*np.exp(-c*(x**d))
def diagPs(x, a,b,c):
    return a+b*np.exp(-c*x)

args = sys.argv

hName = 'hh_'+args[1]+'_'+args[2]+'_'+args[3]+'_'+args[4]+'_'+args[5]+'.dat'
rName = 'hr_'+args[1]+'_'+args[2]+'_'+args[3]+'_'+args[4]+'_'+args[5]+'.dat'
outName = 'pdf_'+args[1]+'_'+args[2]+'_'+args[3]+'_'+args[4]+'_'+args[5]


#load histogram file
df=pd.read_csv(hName,sep=',',header=None,comment='#')
X=df.drop([0,1],axis=1)
X.columns = range(X.shape[1])
width = df[0]
width=np.mean(width) #actually all equal
#average bin -> reduces error on the bin
x=np.mean(X,axis=0)
hgt=(x/np.sum(x))/width
thr=np.arange(X.shape[1])*width
#save results
np.savetxt(outName+'.dat', np.concatenate((thr.reshape(-1,1)[1:], hgt.values.reshape(-1,1)[1:]), axis=1))

#fitting
import scipy.optimize as opt

p,_=opt.curve_fit(Ps, thr, hgt.values)
print('Parameters:',p)

#plot results
plt.scatter(thr, hgt, label='Data', edgecolors='black', c='xkcd:green')
plt.plot(np.linspace(0.01,4), Ps(np.linspace(0.01,4), p[0], p[1],p[2], p[3]), c='xkcd:orangered', label='Fitted density', lw=3)
if(int(args[4])==0):
	plt.plot(np.linspace(0.01,4), Ps(np.linspace(0.01,4), 32/(np.pi**2), 2,4/np.pi,2), c='xkcd:darkblue', label='Wigner Surmise', lw=3)
plt.xlabel('Spacing')
plt.ylabel('P(s)')
plt.legend()
plt.savefig(outName+'.png')
plt.show()



R=pd.read_csv(rName,sep=',',header=None,comment='#')
R.columns = range(R.shape[1])
r = np.mean(R,axis=1)
print(np.mean(r),'+/-',np.std(r))









