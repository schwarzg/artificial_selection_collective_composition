import numpy as np

def community_function(w,m):
	N=w+m
	f=np.divide(m,N)
	return r+s*f


mu=1e-4
r=0.5
s=5e-2
w0=1000
m0=0
tmax=10
nsample=1000
folder="data/time/%d_%d_%s_%s_%s_%s/"%(w0,m0,r,s,mu,tmax)
t=5.0


w,m=np.loadtxt(folder+"/%.1f.dat"%t,unpack=True)
print(w,m)
c=community_function(w,m)

import matplotlib.pyplot as plt

ax=plt.axes((0.1,0.1,0.4,0.3))
ax.hist(c,bins=30)
ax.set_yscale('log')
plt.show()	
