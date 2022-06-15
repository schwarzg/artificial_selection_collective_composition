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
ax.annotate('a',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
ax,_=violinplot_from_histogram(ax,sel_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
ax,_=violinplot_from_histogram(ax,ext_pdf[mask],fbins_t[mask],positions=f0[mask],side='right',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C1',alpha=0.5,marker=4,showmeans=False)
ax.scatter(f0[mask],(sel_median-f0)[mask],c='C0',marker=5)
ax.scatter(f0[mask],(ext_median-f0)[mask],c='C1',marker=4)
ax.hlines(0,0,1,colors='black',ls=':')
ax.set_xlim(xmin=0,xmax=1)
ax.set_ylim(ymin=-0.0401,ymax=0.0601)
ax.set_ylabel(r'$\Delta f^*=f^*_{k+1}-f^*_k$')

#legend plot
lx=plt.axes((0.49,0.857,0.1,0.13))
lxrange=np.arange(-3,3.01,0.03) 
lpx=[st.norm.pdf(lxrange[:-1])]
lx,_=violinplot_from_histogram(lx,lpx,lxrange,color='C0',alpha=0.5,side='left',marker=5) 
lx,_=violinplot_from_histogram(lx,lpx,lxrange,color='C1',alpha=0.5,side='right',marker=4)
#lx.axis('off')
lx.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
lx.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
lx.set_xlim(xmin=-0.4,xmax=0.4)
lx.annotate('Sim',xy=(-0.37,1.9),annotation_clip=False,color='C0') 
lx.annotate('Th',xy=(0.07,1.9),color='C1') 

bx=plt.axes((0.1,0.1,0.5,0.4))
bx.annotate('b',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
bx.plot(f0,sel_median-f0,c='C0',label='Simulation')
bx.plot(f0,ext_median-f0,c='C1',label='Theory')
bx.hlines(0,0,1,colors='black',ls=':')
bx.vlines([fl_ext,fu_ext],-0.05,0.05,colors='black',ls='--')
bx.annotate(r'$f^L$',xy=(fl_ext,-0.0215),xytext=(fl_ext-0.1,-0.019),arrowprops=dict(arrowstyle='->'))
bx.annotate(r'$f^U$',xy=(fu_ext,-0.022),xytext=(fu_ext+0.06,-0.019),arrowprops=dict(arrowstyle='->'))
bx.set_ylim(ymin=-0.022,ymax=0.022)
bx.set_xlim(xmin=0,xmax=1)
bx.set_xlabel(r'Selected Frac. in Cycle $k$, $f^*_k$')
bx.set_ylabel(r'$\mathrm{Median}(\Delta f^*)$')
bx.legend(handlelength=0.6,frameon=False,loc=8)
#plt.savefig('One_step_diagram_g%d_s%s_v2.svg'%(ncomm,s),dpi=300,bbox_inches='tight',format='svg')
plt.show()
