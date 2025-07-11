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
folder2="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,1000)
foldere="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix2"%(N0,r,mu,s,ncomm,1000)

#Simulation
sel_opdf=np.loadtxt(foldere+"_sim.dat")
sel_pdf=np.loadtxt(foldere+"_sel.dat")
sel_median=np.loadtxt(foldere+"_smedian.dat")
sel_mean=np.loadtxt(foldere+"_smean.dat")

#Theory
ext_opdf=np.loadtxt(foldere+"_th.dat")
ext_pdf=np.loadtxt(foldere+"_ext.dat")
ext_median=np.loadtxt(foldere+"_emedian.dat")
ext_mean=np.loadtxt(foldere+"_emean.dat")

fbins_t=[]
fbins_to=[]
fbins_so=[]
fbins_s=[]
f0=np.arange(0.01,0.99,0.01)
f0e=np.arange(0.01,0.99,0.01)
for fsel in f0: 
	#print("Frequency ",fsel)
	fbins_s.append(np.linspace(np.maximum(0,fsel-0.15),np.minimum(1,fsel+0.15),30))
	fbins_so.append(np.linspace(np.maximum(0,fsel-0.15+0),np.minimum(1,fsel+0.15+0.1),30))

for fsel in f0e: 
	fbins_t.append(np.linspace(np.maximum(0,fsel-0.15),np.minimum(1,fsel+0.15),100))
	fbins_to.append(np.linspace(np.maximum(0,fsel-0.15),np.minimum(1,fsel+0.15+0.1),100))

fbins_s=np.array(fbins_s)
fbins_t=np.array(fbins_t)
fbins_to=np.array(fbins_to)
fbins_so=np.array(fbins_so)
mask=range(9,len(f0e),10)
maske=range(9,len(f0e),10)
print(f0e[maske])

#calculate fufl
extmean_interp=itp.interp1d(f0e,ext_median-f0e)
fl_ext=0#opt.root(extmean_interp,0.3).x
fu_ext=0#opt.root(extmean_interp,0.8).x
print(fl_ext,fu_ext)
print(np.shape(sel_pdf[mask]))
##vioiline plot
ax=plt.axes((0.1,0.6,0.5,0.4))
ax.annotate('a',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
#ax.annotate(r'Mutant frequency of the selected Adult $\Psi(f_{k+1}-f_k^*|f_k^*)$',(-0.20,1.05),xycoords='axes fraction')
ax,_=violinplot_from_histogram(ax,sel_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
ax,_=violinplot_from_histogram(ax,ext_pdf[maske],fbins_t[maske],positions=f0e[maske],side='right',width=f0e[maske][1]-f0e[maske][0],yoffset=f0e[maske],color='C1',alpha=0.5,marker=4,showmeans=False)
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


#fig B
bx=plt.axes((0.1,0.05,0.5,0.4))
bx.annotate('b',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')


mask=range(9,len(f0),20)
maske=range(9,len(f0e),20)

bx,_=violinplot_from_histogram(bx,sel_opdf[mask],fbins_so[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=1.,marker=5,showmeans=False,fc='white',ec='C0',lw=2)
bx,_=violinplot_from_histogram(bx,ext_opdf[maske],fbins_to[maske],positions=f0e[maske],side='right',width=f0e[maske][1]-f0e[maske][0],yoffset=f0e[maske],color='C1',alpha=1.,marker=4,showmeans=False,fc='white',ec='C1',lw=2)
bx,_=violinplot_from_histogram(bx,sel_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
bx,_=violinplot_from_histogram(bx,ext_pdf[maske],fbins_t[maske],positions=f0e[maske],side='right',width=f0e[maske][1]-f0e[maske][0],yoffset=f0e[maske],color='C1',alpha=0.5,marker=4,showmeans=False)

bx.hlines(0,0,1,colors='black',ls=':')
bx.set_xlabel(r"Selected F frequency  $f^*_k$")
#bx.set_ylabel(r"$\bar{f}_{k+1}-f^*_k$")
#bx.set_ylabel(r"$\sigma_{{f}_{k+1}}$")
bx.set_ylabel(r'$f_{k+1}^*-f^*_k$')
bx.legend(frameon=False)
#bx.set_xlim(xmin=0,xmax=1)
bx.set_ylim(ymin=-0.0401,ymax=0.1601)



formatter='svg'
plt.savefig('figures/FigS5_fix.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
