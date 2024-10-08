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
nens= 300 # 1000 for +, 300 for -
ncycle=10
tcycle=np.log(ncomm+1)/r

from one_step_functions import *
from custom_plot import *

##Load data
folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)

#Simulation
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
	#print("Frequency ",fsel)
	fbins_s.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))
	fbins_t.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))

fbins_s=np.array(fbins_s)
fbins_t=np.array(fbins_t)
mask=range(10,len(f0),10)

#calculate fufl
extmean_interp=itp.interp1d(f0,ext_median-f0)
fl_ext=opt.root(extmean_interp,0.3).x
fu_ext=opt.root(extmean_interp,0.8).x
print(fl_ext,fu_ext)

##vioiline plot
ax=plt.axes((0.1,0.6,0.5,0.4))
#ax.annotate('a',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
#ax.annotate(r'Mutant frequency of the selected Adult $\Psi(f_{k+1}-f_k^*|f_k^*)$',(-0.20,1.05),xycoords='axes fraction')
ax,_=violinplot_from_histogram(ax,sel_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
ax,_=violinplot_from_histogram(ax,ext_pdf[mask],fbins_t[mask],positions=f0[mask],side='right',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C1',alpha=0.5,marker=4,showmeans=False)
#ax.scatter(f0[mask],(sel_median-f0)[mask],c='C0',marker=5)
#ax.scatter(f0[mask],(ext_median-f0)[mask],c='C1',marker=4)
ax.hlines(0,0,1,colors='black',ls=':')
ax.set_xlim(xmin=0,xmax=1)
ax.set_ylim(ymin=-0.0401,ymax=0.0601)
ax.set_ylabel(r'$f_{k+1}^*-f^*_k$')
ax.set_xlabel(r'Selected F frequency in cycle $k$, $f^*_k$')

#legend plot
lx=plt.axes((0.49,0.857,0.1,0.13))
lxrange=np.arange(-3,3.01,0.03) 
lpx=[st.norm.pdf(lxrange[:-1])]
lx,_=violinplot_from_histogram(lx,lpx,lxrange,color='C0',alpha=0.5,side='left',showmeans=False) 
lx,_=violinplot_from_histogram(lx,lpx,lxrange,color='C1',alpha=0.5,side='right',showmeans=False)
#lx.axis('off')
lx.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
lx.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
lx.set_xlim(xmin=-0.4,xmax=0.4)
lx.annotate('Sim',xy=(-0.37,1.9),annotation_clip=False,color='C0') 
lx.annotate('Th',xy=(0.07,1.9),color='C1') 
'''
#fig B
bx=plt.axes((0.72,0.6,0.5,0.4))
bx.annotate('b',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')

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
'''



formatter='svg'
plt.savefig('figures/FigS5.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
