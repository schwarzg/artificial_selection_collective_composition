import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.interpolate as itp
import scipy.optimize as opt
import sys
sys.setrecursionlimit(10000)
#parameter prepare
mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
rhat=0
nens=100
ncycle=10
tcycle=np.log(ncomm+1)/r


from one_step_functions import *
from selection import select
from selection import community_function as cf
from selection import score_function as sf
from reproduction import hypergeometric_reproduce
from scipy.signal import savgol_filter

#Start from give initial state
Nf=(ncomm+2)*N0		#Nf is set to be (ncomm+1)*N0
fbins=np.linspace(0,1,100)

#frequency loop
fls=np.array([])
fle=np.array([])
fus=np.array([])
fue=np.array([])
nens=1000
for ncomm in [4,6,8,10,12,14,16,18,20,22]:
	folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
	seldat=np.loadtxt(folder+"_smedian.dat")
	extdat=np.loadtxt(folder+"_emedian.dat")
	f0=np.arange(0.00,0.99,0.01)
	
	hatsel=savgol_filter(seldat-f0,11,2)	#smoothing
	selmean_interp=itp.interp1d(f0,hatsel)
	extmean_interp=itp.interp1d(f0,extdat-f0)
	fl_sel=opt.root(selmean_interp,0.3).x
	try:
		fu_sel=opt.root(selmean_interp,0.8).x
	except:
		fu_sel=opt.root(selmean_interp,0.6).x
	fl_ext=opt.root(extmean_interp,0.3).x
	fu_ext=opt.root(extmean_interp,0.8).x
	print(ncomm,fl_sel,fu_sel,fl_ext,fu_ext)
	fls=np.append(fls,fl_sel)
	fus=np.append(fus,fu_sel)
	fle=np.append(fle,fl_ext)
	fue=np.append(fue,fu_ext)

#N0 data
ncomm=10
fls2=np.array([])
fle2=np.array([])
fus2=np.array([])
fue2=np.array([])
N0s=[800,900,1000,2000,4000,6000]
for N0 in N0s:
	folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
	seldat=np.loadtxt(folder+"_smedian.dat")
	extdat=np.loadtxt(folder+"_emedian.dat")
	f0=np.arange(0.00,0.99,0.01)[:len(seldat)]
	
	hatsel=savgol_filter(seldat-f0,11,2)	#smoothing
	selmean_interp=itp.interp1d(f0,hatsel,fill_value='extrapolate')
	extmean_interp=itp.interp1d(f0,extdat-f0,fill_value='extrapolate')
	
	#search range when fl is approachable
	left1=0.01
	right1=0.02
	while selmean_interp(left1)*selmean_interp(right1)>0 and right1<1.0:	 #until the sign is different
		right1=right1+0.01
	left2=0.51
	right2=0.52
	while selmean_interp(left2)*selmean_interp(right2)>0 and right2<1.0:	 #until the sign is different
		right2=right2+0.01
	print(right1,right2)
	if np.abs(right1-1)<1e-8:
		fl_sel=0.5
		fu_sel=0.5
	else:
		fl_sel=opt.root(selmean_interp,right1).x
		fu_sel=opt.root(selmean_interp,right2).x

	left1=0.1
	right1=0.02
	while extmean_interp(left1)*extmean_interp(right1)>0 and right1<1.0:	 #until the sign is different
		right1=right1+0.01
	left2=0.51
	right2=0.52
	while extmean_interp(left2)*extmean_interp(right2)>0 and right2<1.0:	 #until the sign is different
		right2=right2+0.01
	print(right1,right2)
	if np.abs(right1-1)<0.1:
		fl_ext=0.5
		fu_ext=0.5
	else:
		fl_ext=opt.root(extmean_interp,right1).x
		fu_ext=opt.root(extmean_interp,right2).x
	print(N0,fl_sel,fu_sel,fl_ext,fu_ext)
	fls2=np.append(fls2,fl_sel)
	fus2=np.append(fus2,fu_sel)
	fle2=np.append(fle2,fl_ext)
	fue2=np.append(fue2,fu_ext)

shadeup=1.1*np.ones(len(N0s))
shadedw=np.zeros(len(N0s))
print(shadeup,shadedw,fue2)
ax=plt.axes((0.1,0.1,0.3,0.6))
ax.annotate('a',(-0.22,1.08),xycoords='axes fraction',fontweight='bold')
ax.set_ylim(ymin=0,ymax=1)
ax.plot(N0s,fus2,c='C0',label='Sim,$f^U$',marker='v',ms=4)
ax.plot(N0s,fue2,c='C1',label='Th,$f^U$',marker='^',ms=4)
ax.fill_between(N0s,shadeup,fue2,color='lightgray')
ax.plot(N0s,fls2,c='C0',label='Sim,$f^L$',marker='v',ms=4,ls='--')
ax.plot(N0s,fle2,c='C1',label='Th,$f^L$',marker='^',ms=4,ls='--')
ax.fill_between(N0s,fle2,shadedw,color='lightgray')
ax.set_xscale('log')
ax.set_xlim(xmin=500,xmax=6000)
ax.fill_between([1,900],[1.1,1.1],color='lightgray')
ax.set_ylabel(r'Mutant Frequency $f$')
ax.set_xlabel(r'Newborn collective size $N_0$')
ax.set_xticks([500,1000,5000])
ax.set_xticklabels([500,1000,5000])
ax.legend(frameon=False)

bx=plt.axes((0.5,0.1,0.3,0.6))
fs=[4,6,8,10,12,14,16,18,20,22]
shadeup=1.1*np.ones(len(fs))
shadedw=np.zeros(len(fs))
bx.annotate('b',(-0.22,1.08),xycoords='axes fraction',fontweight='bold')
bx.set_ylim(ymin=0,ymax=1)
bx.plot(fs,fus,c='C0',label='Sim,$f^U$',marker='v',ms=4)
bx.plot(fs,fue,c='C1',label='Th,$f^U$',marker='^',ms=4)
bx.fill_between(fs,shadeup,fue,color='lightgray')
bx.plot(fs,fls,c='C0',label='Sim,$f^L$',marker='v',ms=4,ls='--')
bx.plot(fs,fle,c='C1',label='Th,$f^L$',marker='^',ms=4,ls='--')
bx.fill_between(fs,fle,shadedw,color='lightgray')
bx.set_xlabel(r'Number of collectives $g$')
bx.set_xlim(4,22)

plt.savefig('Fig5.svg',dpi=300,bbox_inches='tight',format='svg')
plt.show()

