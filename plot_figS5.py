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
s=-3e-2
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
f0=np.arange(0.01,0.99,0.01)
for fsel in f0: 
	#print("Frequency ",fsel)
	fbins_s.append(np.linspace(np.maximum(0,fsel-0.15),np.minimum(1,fsel+0.15),100))
	fbins_t.append(np.linspace(np.maximum(0,fsel-0.15),np.minimum(1,fsel+0.15),100))

fbins_s=np.array(fbins_s)
fbins_t=np.array(fbins_t)
mask=range(9,len(f0),10)

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
ax.set_ylabel(r'$f-f^*_k$')

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
bx.set_ylabel(r'$\mathrm{Median}[\Psi(f-f^*_k|f^*_k)]$')
bx.legend(handlelength=0.6,frameon=False,loc=8)


#heat map diagram
ncycle=1000
dat=np.loadtxt("data/ens/N0%s_r%s_s%s_mu%s_g%s_ncycle%d_diagram_fig2.txt"%(N0,r,s,mu,ncomm,ncycle))

rhats=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
mbars=np.array([0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])/N0

cx=plt.axes((0.75,0.3,0.5,0.5))
cx.annotate('c',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
frange=np.arange(0.1,1.0,0.1)
heatmap=cx.pcolormesh(rhats,mbars,dat,cmap='Blues_r',shading='flat')
heatmap.set_clim(0,0.05)
cbar=plt.colorbar(heatmap,extend='max')
cbar.set_label(r'$\langle| f^*_{k_\mathrm{final}}-\hat{f}|\rangle/\hat{f}$')
cx.set_xlim(xmin=0,xmax=1)
cx.set_ylim(ymin=0,ymax=1)
cx.set_xlabel(r'Initial Frequency $\bar{f}_o$')
cx.set_ylabel(r'Target Frequency $\hat{f}$')

'''
#Suppoting information
cx.scatter([0.05],[0.85],c='C1',marker='^')
cx.scatter([0.05],[0.50],c='C2',marker='^')
cx.scatter([0.05],[0.15],c='C3',marker='^')
cx.scatter([0.5],[0.85],c='C1',marker='v')
cx.scatter([0.5],[0.50],c='C2',marker='v')
cx.scatter([0.5],[0.15],c='C3',marker='v')
fl=0.28580445
fu=0.68687363
cx.hlines(fl,0,fl,colors='black',ls='--')
cx.vlines(fl,0,fl,colors='black',ls='--')
cx.hlines(fu,0,1,colors='black',ls='--')
cx.annotate(r'$f^L$',xy=(0,fl),xytext=(0.10,fl+0.10),arrowprops=dict(arrowstyle='->'))
cx.annotate(r'$f^L$',xy=(fl,0),xytext=(fl+0.10,0.10),arrowprops=dict(arrowstyle='->'))
cx.annotate(r'$f^U$',xy=(0,fu),xytext=(0.10,fu-0.10),arrowprops=dict(arrowstyle='->'))

'''
############################################
# for negative s
############################################

print(fl_ext,fu_ext)
fl=fl_ext#0.31601399
fu=fu_ext#0.71544258
cx.hlines(fu,1,fu,colors='black',ls='--')
cx.vlines(fu,1,fu,colors='black',ls='--')
cx.hlines(fl,0,1,colors='black',ls='--')
cx.annotate(r'$f^U$',xy=(1,fu),xytext=(1-0.13,fu-0.13),arrowprops=dict(arrowstyle='->'))
cx.annotate(r'$f^U$',xy=(fu,1),xytext=(fu-0.13,1-0.13),arrowprops=dict(arrowstyle='->'))
cx.annotate(r'$f^L$',xy=(0,fl),xytext=(0.10,fl+0.10),arrowprops=dict(arrowstyle='->'))
#cx.set_xticks([0.0,0.2,fl,0.4,0.6,0.8,1.0])
#ax.set_xticklabels(['0.0','0.2',r'$f^L$','0.4','0.6','0.8','1.0'])

formatter='svg'
plt.savefig('FigS5.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
