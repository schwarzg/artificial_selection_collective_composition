import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "arial"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.rcParams["mathtext.fontset"]='stix'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'STIXGeneral:italic:bold'
plt.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'


from analytic_results import *

mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
rhat=0.05
nens=300
ncycle=10
tcycle=np.log(ncomm+1)/r

def get_f(w,m):
	return np.divide(m,w+m)

barN02=[]
sigw2=[]
sigm2=[]
sigf2=[]
N0s= np.arange(500,5000,100)[::-1]
for N0 in N0s:
	m0=int(N0*0.3)#np.clip(np.random.poisson(200,size=nens),0,1000)
	w0=N0-m0
	c0=get_f(w0,m0)
	barw=barw_th(np.mean(w0),tcycle,r,mu)	
	sig2w=sig2w_th(np.mean(w0),np.var(w0),tcycle,r,mu)
	barm=barm_th(np.mean(m0),np.mean(w0),tcycle,r,mu,s)
	sig2m=sig2m_th(np.mean(m0),np.var(m0),np.mean(w0),np.var(w0),-np.var(m0),tcycle,r,mu,s)	
	sig2c=sig2c_th_v2(np.mean(c0),np.var(c0),N0,tcycle,r,mu,s)

	barN02.append((barw+barm)**2)
	sigw2.append(sig2w)
	sigm2.append(sig2m)
	sigf2.append(sig2c)

#ax=plt.axes((0.1,0.1,0.3,0.3))
#bx=plt.axes((0.53,0.1,0.3,0.3))
'''
ax.plot(N0s,sigw2,label=r'$\sigma_w^2(\tau)$')	
ax.plot(N0s,sigm2,label=r'$\sigma_m^2(\tau)$')
ax.plot([1500,3000],[45000,90000],ls=':',c='black')
ax.annotate(r'$\sim N_0$',(2000,45000))
ax.set_xlabel(r'Newborn collective size $N_0$')	
ax.set_ylabel(r'Variances')	
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(frameon=False,handlelength=0.7)
ax.set_xticks([500,1000,5000])
ax.set_xticklabels([500,1000,5000])
bx.plot(N0s,barN02)	
bx.plot([1500,3000],[150000000,600000000],ls=':',c='black')
bx.annotate(r'$\sim N_0^2$',(2000,150000000))
bx.set_xscale('log')
bx.set_yscale('log')
bx.set_xlabel(r'Newborn collective size $N_0$')	
bx.set_ylabel(r'$\overline{N}(\tau)^2$')	
bx.set_xticks([500,1000,5000])
bx.set_xticklabels([500,1000,5000])
'''
bx=plt.axes((0.1,0.1,0.3,0.3))
bx.annotate('a',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')

from analytic_results import barc_th_v2,sig2c_th_v2

f0s = np.arange(0,1+0.01*0.5,0.01)
barcs=[]
sig2cs=[]
for f0 in f0s:
    barc=barc_th_v2(f0,N0,tcycle,r,mu,s)
    barcs.append(barc)
    sig2c=sig2c_th_v2(f0,f0*(1-f0)/N0,N0,tcycle,r,mu,s)
    sig2cs.append(sig2c)

bx.plot(f0s,barcs-f0s,label=r"$\bar{f}_{k+1}-f^*_k$")
bx.plot(f0s,np.sqrt(sig2cs),label=r"$\sigma_{{f}_{k+1}}$")

bx.set_xlabel(r"Selected F frequency  $f^*_k$")
#bx.set_ylabel(r"$\bar{f}_{k+1}-f^*_k$")
#bx.set_ylabel(r"$\sigma_{{f}_{k+1}}$")
bx.set_ylabel(r"Values in frequency $f$")
bx.legend(frameon=False)

cx=plt.axes((0.55,0.1,0.3,0.3))
cx.annotate('b',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')
cx.plot(1/N0s,sigf2)	
cx.plot([1/1000,1/2000],[0.00016,0.00008],ls=':',c='black')
cx.annotate(r'$\sim 1/N_0$',(1/1000,0.00012))
cx.set_xscale('log')
cx.set_yscale('log')
cx.set_xlabel(r'Newborn collective size $1/N_0$')	
cx.set_ylabel(r'$\sigma_f^2(\tau)$')	
cx.set_xticks([1/1000,1/5000])
cx.set_xticklabels(['1/1000','1/5000'])
formatter='svg'
plt.savefig('figures/FigS6.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
