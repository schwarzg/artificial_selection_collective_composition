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
	#print("Frequency ",fsel)
	fbins_s.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))
	fbins_t.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))

fbins_s=np.array(fbins_s)
fbins_t=np.array(fbins_t)
mask=[9,49,89]#=range(10,len(f0),10)

#calculate fufl
extmean_interp=itp.interp1d(f0,ext_median-f0)
fl_ext=opt.root(extmean_interp,0.3).x
fu_ext=opt.root(extmean_interp,0.8).x
print(fl_ext,fu_ext)

##vioiline plot
ax=plt.axes((0.1,0.85,0.5,0.35))
ax.annotate('a',(-0.25,1.07),xycoords='axes fraction',fontweight='bold')
ax.annotate(r'Median of Mutant frequency of the selected Adult, $\Psi(f^*_{k+1}=f|f_k^*)$',(-0.20,1.07),xycoords='axes fraction')
#bx,_=violinplot_from_histogram(bx,ext_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
#bx,_=violinplot_from_histogram(ax,ext_pdf[mask],fbins_t[mask],positions=f0[mask],side='right',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C1',alpha=0.5,marker=4,showmeans=False)
#ax.scatter(f0[mask],(sel_median-f0)[mask],c='C0',marker=5)
#ax.scatter(f0[mask],(ext_median-f0)[mask],c='C1',marker=4)
ax.hlines(0,0,1,colors='black',ls=':')
ax.set_xlim(xmin=0,xmax=1)
#ax.set_ylim(ymin=-0.0401,ymax=0.0601)
ax.set_ylim(ymin=-0.022,ymax=0.022)
ax.set_ylabel(r'$f-f^*_k$')

'''
#legend plot
lx=plt.axes((0.492,1.092,0.1,0.11))
lxrange=np.arange(-3,3.01,0.03) 
lpx=[st.norm.pdf(lxrange[:-1])]
lx,_=violinplot_from_histogram(lx,lpx,lxrange,color='C0',alpha=0.5,side='left',showmeans=False) 
lx,_=violinplot_from_histogram(lx,lpx,lxrange,color='C1',alpha=0.5,side='right',showmeans=False)
#lx.axis('off')
lx.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
lx.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
lx.set_xlim(xmin=-0.4,xmax=0.4)
lx.annotate('Sim',xy=(-0.37,1.82),annotation_clip=False,color='C0') 
lx.annotate('Th',xy=(0.07,1.82),color='C1') 
'''

#median plot
bx=plt.axes((0.1,0.38,0.5,0.35))
#bx.annotate('b',(-0.25,1.07),xycoords='axes fraction',fontweight='bold')
#bx.annotate(r'Median of $\Psi$',(-0.20,1.07),xycoords='axes fraction')
ax.plot(f0,sel_median-f0,c='C0',label='Simulation')
ax.plot(f0,ext_median-f0,c='C1',label='Theory')
#bx.hlines(0,0,1,colors='black',ls=':')
ax.vlines([fl_ext,fu_ext],-0.05,0.05,colors='black',ls='--')
ax.annotate(r'$f^L$',xy=(fl_ext,-0.0215),xytext=(fl_ext-0.1,-0.019),arrowprops=dict(arrowstyle='->'))
ax.annotate(r'$f^U$',xy=(fu_ext,-0.022),xytext=(fu_ext+0.06,-0.019),arrowprops=dict(arrowstyle='->'))
#bx.set_ylim(ymin=-0.022,ymax=0.022)
#bx.set_xlim(xmin=0,xmax=1)
#bx.set_xlabel(r" Selected mutant frequency" "\n" "of the selected Adult in cycle $k$, $f^*_k$")
#bx.set_ylabel(r'$f-f^*_k$')
#bx.legend(handlelength=0.6,frameon=False,loc=8)

#Difference plot
#cx=plt.axes((0.1,0,0.5,0.25))
#dx=plt.axes((0.4,0,0.2,0.25))
bx.annotate('b',(-0.25,1.07),xycoords='axes fraction',fontweight='bold')
bx.annotate(r'Distribution of frequency difference in simulation',(-0.20,1.07),xycoords='axes fraction')
bx.set_ylabel(r'Frequency difference')


idx1=9
idx2=49
idx3=89
bx.hlines(0,0,1,colors='black',ls=':')
c=['C0','C1','C2']
boxes=bx.boxplot([sim_raw[idx1]-f0[idx1],sel_raw[idx1]-np.mean(sim_raw[idx1]),sel_raw[idx1]-f0[idx1]],positions=[0.04,0.1,0.16],widths=0.05,sym='')
for box,med,col in zip(boxes['boxes'],boxes['medians'],c):
	box.set_color(col)
	med.set_color(col)
boxes=bx.boxplot([sim_raw[idx2]-f0[idx2],sel_raw[idx2]-np.mean(sim_raw[idx2]),sel_raw[idx2]-f0[idx2]],positions=[0.44,0.5,0.56],widths=0.05,sym='')
for box,med,col in zip(boxes['boxes'],boxes['medians'],c):
	box.set_color(col)
	med.set_color(col)
boxes=bx.boxplot([sim_raw[idx3]-f0[idx3],sel_raw[idx3]-np.mean(sim_raw[idx3]),sel_raw[idx3]-f0[idx3]],positions=[0.84,0.9,0.96],widths=0.05,sym='')
for box,med,col in zip(boxes['boxes'],boxes['medians'],c):
	box.set_color(col)
	med.set_color(col)
#bx,_=violinplot_from_histogram(bx,ext_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
bx.set_xlim(xmin=0,xmax=1)
bx.set_xticks([0.1,0.5,0.9])
bx.set_xticklabels(['0.1','0.5','0.9'])
bx=statistical_significant(bx,sim_raw[idx1]-f0[idx1],sim_raw[idx2]-f0[idx2],[f0[idx1]-0.05,f0[idx2]-0.05],upside=False)
bx=statistical_significant(bx,sel_raw[idx1]-np.mean(sim_raw[idx1]),sel_raw[idx2]-np.mean(sim_raw[idx2]),[f0[idx1]+0.01,f0[idx2]+0.01],upside=False,offset=0.1)


bx.annotate('Intra-collective',(0.03,0.80),xycoords='axes fraction',rotation=15,c=c[0])
bx.annotate('Inter-collective',(0.08,0.70),xycoords='axes fraction',rotation=15,c=c[1])
bx.annotate('Total',(0.14,0.61),xycoords='axes fraction',rotation=15,c=c[2])
bx.set_xlabel(r" Selected mutant frequency" "\n" "of the selected Adult in cycle $k$, $f^*_k$")

#legend for boxplot
#lcx=plt.axes((0.48,0.1,0.1,0.13,))
#lcx.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
#lcx.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
#plt.savefig('Fig4_v2.svg',dpi=300,bbox_inches='tight',format='svg')
#plt.savefig('Fig4_v2.png',dpi=300,bbox_inches='tight')
plt.show()
