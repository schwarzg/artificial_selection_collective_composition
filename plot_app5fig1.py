import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'DejaVu Sans:italic:bold'
plt.rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
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
nens= 1000 # 1000 for +, 300 for -
ncycle=10
tcycle=np.log(ncomm+1)/r

from one_step_functions import *
from custom_plot import *
##Load data
folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)

#Simulation
sim_raw=np.loadtxt(folder+"_simraw.dat")
sel_raw=np.loadtxt(folder+"_selraw.dat")
sel_pdf=np.loadtxt(folder+"_sel.dat")
sel_median=np.loadtxt(folder+"_smedian.dat")
sel_mean=np.loadtxt(folder+"_smean.dat")

#Theory
ext_pdf=np.loadtxt(folder+"_ext.dat")
ext_median=np.loadtxt(folder+"_emedian.dat")
ext_mean=np.loadtxt(folder+"_emean.dat")

fbins_t=[]
fbins_s=[]
f0=np.arange(0.00,0.99,0.01)
for fsel in f0: 
	fbins_s.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))
	fbins_t.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))

fbins_s=np.array(fbins_s)
fbins_t=np.array(fbins_t)
mask=[9,49,89]#=range(10,len(f0),10)

#calculate fufl
extmean_interp=itp.interp1d(f0,ext_median-f0)
fl_ext=opt.root(extmean_interp,0.3).x
fu_ext=opt.root(extmean_interp,0.8).x
import scipy as sc
import scipy.optimize as opt
import scipy.special as spc
import scipy.stats as st

#define D function
def get_Dfunc(zeta,r,s,mu,ncomm,N0):
	varf0=zeta*(1-zeta)/N0
	barc=barc_th_v2(zeta,N0,tcycle,r,mu,s)
	sig2c=sig2c_th_v2(zeta,varf0,N0,tcycle,r,mu,s)
	unitnorm=st.norm()
	phi1=unitnorm.ppf(float(np.log(2)/ncomm))
	phi2=unitnorm.ppf(float(1/(ncomm*np.exp(1))))
	#D=barc+np.sqrt(sig2c)*(phi1-np.log(np.log(2))*(phi2-phi1))-zeta
	D=barc+np.sqrt(sig2c)*(phi1)-zeta
	return D

zetas=np.arange(0,1,0.01)
ss=np.arange(0.02,0.060005,0.00001)
fls=[]
fus=[]
lastloss=None
for s in ss:
    Ds=get_Dfunc(zetas,r,s,mu,ncomm,N0)
    extmean_interp=itp.interp1d(zetas,Ds)
    #plt.plot(zetas,Ds)
    #plt.hlines(0,1,0)
    #plt.show()
    
    solved=False
    lsolve=False
    usolve=False
    sl=0.3
    su=0.7
    fl_th=None
    fu_th=None
    while not solved:
        if not np.sum(Ds[:1]*Ds[1:]<0): # no passing 
            fls.append(None)
            fus.append(None)
            solved=True
            print(s,' is solutionless')
            lastloss=s
            continue
        if not lsolve:
            try:
                fl_th=opt.root(extmean_interp,sl).x
            except:
                sl=max(sl-0.05,0)
            else:
                fls.append(fl_th[0])
                lsolve=True

        if not usolve:
            try:
                fu_th=opt.root(extmean_interp,su).x
            except:
                su=min(su+0.05,1)
            else:
                fus.append(fu_th[0])
                usolve=True
        solved=usolve and lsolve

print(fls,fus)
Ds=get_Dfunc(zetas,r,lastloss,mu,ncomm,N0)

plt.plot(zetas,Ds+0.00001)
plt.hlines(0,1,0)
plt.show()

args=np.argwhere(fls)[:,0].astype(int)
print(args)

ax=plt.axes((0.1,0.1,0.4,0.3))
ax.set_xlabel(r'Selecive advantage $\omega$')
ax.set_ylabel(r'Mutant freuqency $f$')
ax.plot(ss,fls,label=r'$f^L$',c='C2',ls='--')
ax.plot(ss,fus,label=r'$f^H$',c='C2')
ax.fill_between(ss[296:],fus[296:],fls[296:],color='gray')
#ax.legend(frameon=False)
ax.set_xlim(xmin=0.01)

#Ds0=get_Dfunc(zetas,r,s,0,ncomm,N0)

#plot position assign

formatter='svg'
#plt.savefig('figures/FigSxx.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
