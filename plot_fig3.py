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
folder2="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix"%(N0,r,mu,s,ncomm,nens)

#Simulation
sim_raw=np.loadtxt(folder+"_simraw.dat")
sel_raw=np.loadtxt(folder+"_selraw.dat")
sel_pdf=np.loadtxt(folder+"_sel.dat")
sel_median=np.loadtxt(folder+"_smedian.dat")
sel_mean=np.loadtxt(folder+"_smean.dat")

#Theory
ext_pdf=np.loadtxt(folder2+"_ext.dat")
ext_median=np.loadtxt(folder2+"_emedian.dat")
ext_mean=np.loadtxt(folder2+"_emean.dat")

fbins_t=[]
fbins_s=[]
f0=np.arange(0.00,0.99,0.01)
for i,fsel in enumerate(f0): 
	fbins_s.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))

f0t=np.arange(0.01,0.99,0.01)
for i,fsel in enumerate(f0t): 
	fbins_t.append(np.linspace(np.maximum(0,fsel-0.05),np.minimum(1,fsel+0.05),30))

fbins_s=np.array(fbins_s)
fbins_t=np.array(fbins_t)
mask=[9,49,89]#=range(10,len(f0),10)

#calculate fufl
extmean_interp=itp.interp1d(f0t,ext_median-f0t)
fl_ext=opt.root(extmean_interp,0.3).x
fu_ext=opt.root(extmean_interp,0.8).x
print(fl_ext,fu_ext)
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
Ds=get_Dfunc(zetas,r,s,mu,ncomm,N0)
#Ds0=get_Dfunc(zetas,r,s,0,ncomm,N0)

#plot position assign
axx=0.1
axy=0.95
bxx=0.1
bxy=0.38
cxx=0.0
cxy=0.0


##median plot
ax=plt.axes((axx,axy,0.5,0.30))
#ax.annotate('a',(-0.25,1.19),xycoords='axes fraction',fontweight='bold')
#ax.annotate(r'Median of selected mutant frequency distribution $\Psi(f^*_{k+1}|f_k^*)$'+'\n'+r'with shift $-f^*_k$',(-0.20,1.07),xycoords='axes fraction')
ax.hlines(0,0,1,colors='black',ls=':')
ax.set_xlim(xmin=0,xmax=1)
ax.set_ylim(ymin=-0.018,ymax=0.018)
ax.set_ylabel(r'$\mathbf{f^*_{k+1}-f^*_k}$')
ax.set_xlabel(r"Selected mutant frequency in cycle $\mathbf{k}$, $\mathbf{f^*_k}$")

ax.plot(f0,(sel_median-f0),c='C0',label='Simulation')
ax.plot(f0t,(ext_median-f0t),c='C1',label='Equation[1]')
ax.plot(zetas,Ds,label=r'$Equation[2]$',c='C2')
#ax.plot(zetas,Ds0,label=r'$D(\zeta)$',c='C3')
ax.legend(frameon=False)
#ax.vlines([fl_ext,fu_ext],-0.05,0.05,colors='black',ls='--')
#ax.annotate(r'$\mathbf{f^L}$',xy=(fl_ext,-0.0215),xytext=(fl_ext-0.1,-0.019),arrowprops=dict(arrowstyle='->'))
#ax.annotate(r'$\mathbf{f^U}$',xy=(fu_ext,-0.022),xytext=(fu_ext+0.06,-0.019),arrowprops=dict(arrowstyle='->'))

#distribution box plot
bx=plt.axes((0.1,0.38,0.5,0.35))
bx.annotate('b',(-0.25,1.07),xycoords='axes fraction',fontweight='bold')
bx.annotate(r'Distribution of frequency changes in each steps',(-0.20,1.07),xycoords='axes fraction')
bx.set_ylabel(r'Frequency difference')


idx1=9
idx2=49
idx3=89
bx.hlines(0,0,1,colors='black',ls=':')
c=['#17becf','#9467bd','black']
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
bx.set_xlim(xmin=0,xmax=1)
bx.set_xticks([0.1,0.5,0.9])
bx.set_xticklabels(['0.1','0.5','0.9'])
bx=statistical_significant(bx,sim_raw[idx1]-f0[idx1],sim_raw[idx2]-f0[idx2],[f0[idx1]-0.05,f0[idx2]-0.05],upside=False)
bx=statistical_significant(bx,sel_raw[idx1]-np.mean(sim_raw[idx1]),sel_raw[idx2]-np.mean(sim_raw[idx2]),[f0[idx1]+0.01,f0[idx2]+0.01],upside=False,offset=0.1)


bx.annotate('Intra-collective',(0.03,0.78),xycoords='axes fraction',rotation=15,c=c[0])
bx.annotate('Inter-collective',(0.07,0.68),xycoords='axes fraction',rotation=15,c=c[1])
bx.annotate('Total',(0.14,0.60),xycoords='axes fraction',rotation=15,c=c[2])
#bx.annotate(r'$f_{k+1,\tau}-f^*_{k}$',(0.03,0.78),xycoords='axes fraction',rotation=15,c=c[0])
#bx.annotate(r'$f^*_{k+1}-\overline{f_{k+1,\tau}}$',(0.03,0.78),xycoords='axes fraction',rotation=15,c=c[1])
bx.set_xlabel(r" Selected mutant frequency in cycle $\mathbf{k}$, $\mathbf{f^*_k}$")

'''
#Conclutsion for Two strain case
#large frequency
cx1=plt.axes((cxx,cxy,0.4,0.25))
cx1.annotate('c',(-0.05,0.8),xycoords='axes fraction',fontweight='bold')
cx1.annotate(r'Success or Fail according to the target composition',(0.00,0.8),xycoords='axes fraction')

#frame
cx1.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
cx1.fill_between([0,1],0.02,0,color='green',alpha=0.3)
cx1.annotate(r'$\mathbf{0}$',xy=(0,-0.02))
cx1.annotate(r'$\mathbf{1}$',xy=(0.95,-0.02))

cx1.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.9
ini1=0.2
h1=0.01
ini2=0.5
h2=0.01
ini3=0.8
h3=0.01
cx1.vlines(tgt,0,0.02,colors='black',ls=':',label='Target') #target
cx1.scatter(ini1,h1,c='black',marker='o',label='Initial') #initial 1
cx1.scatter(ini2,h2,c='black',marker='o') #initial 2
cx1.scatter(ini3,h3,c='black',marker='o') #initial 2

cx1.annotate('',xy=(ini2,h2),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
cx1.annotate('',xy=(ini3,h3),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
cx1.annotate('',xy=(tgt,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))

cx1.vlines([0.3,0.7],0.005,-0.005,colors='black')
cx1.annotate(r'$\mathbf{f^L}$',xy=(0.28,-0.02))
cx1.annotate(r'$\mathbf{f^U}$',xy=(0.68,-0.02))

cx1.set_xlim(xmin=0,xmax=1)

cx1.axis('off')

#mid frequency
cx2=plt.axes((cxx,cxy-0.15,0.4,0.25))

#frame
cx2.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
cx2.fill_between([0,1],0.02,0,color='red',alpha=0.3)
cx2.annotate(r'$\mathbf{0}$',xy=(0,-0.02))
cx2.annotate(r'$\mathbf{1}$',xy=(0.95,-0.02))

cx2.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.5
ini1=0.2
h1=0.010
ini2=0.5
h2=0.01
ini3=0.8
h3=0.01
cx2.vlines(tgt,0,0.02,colors='black',ls=':',label='Target') #target
cx2.scatter(ini1,h1,c='black',marker='o',label='Initial(low)') #initial 1
cx2.scatter(ini2,h2,c='black',marker='o',label='Initial(mid)') #initial 2
cx2.scatter(ini3,h3,c='black',marker='o',label='Initial(high)') #initial 2

cx2.annotate('',xy=(ini2,h2),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
cx2.annotate('',xy=(0.7,h2),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
cx2.annotate('',xy=(0.7,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))

cx2.vlines([0.3,0.7],0.005,-0.005,colors='black')
cx2.annotate(r'$\mathbf{f^L}$',xy=(0.28,-0.02))
cx2.annotate(r'$\mathbf{f^U}$',xy=(0.68,-0.02))

cx2.set_xlim(xmin=0,xmax=1)
cx2.axis('off')

#small frequency
cx3=plt.axes((cxx,cxy-0.3,0.4,0.25))

#frame
cx3.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
cx3.fill_between([0,0.3],0.025,0,color='green',alpha=0.3, label='success')
cx3.fill_between([0.3,1],0.025,0,color='red',alpha=0.3,label='fail')
cx3.annotate(r'$\mathbf{0}$',xy=(0,-0.02))
cx3.annotate(r'$\mathbf{1}$',xy=(0.95,-0.02))
cx3.annotate('Mutant frequency',xy=(0.5,0.2),xycoords='axes fraction',ha='center',va='top')

cx3.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.05
ini1=0.2
h1=0.015
ini2=0.5
h2=0.015
ini3=0.8
h3=0.015

cx3.scatter(ini1,h1,c='black',marker='o',label='Initial') #initial 1
cx3.scatter(ini2,h2,c='black',marker='o') #initial 2
cx3.scatter(ini3,h3,c='black',marker='o') #initial 2
cx3.vlines(tgt,-0.005,0.025,colors='black',ls=':',label='Target') #target

cx3.annotate('',xy=(tgt,h1),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
cx3.annotate('',xy=(0.7,h2),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
cx3.annotate('',xy=(0.7,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))

cx3.vlines([0.3,0.7],0.005,-0.005,colors='black')
cx3.annotate(r'$\mathbf{f^L}$',xy=(0.28,-0.02))
cx3.annotate(r'$\mathbf{f^U}$',xy=(0.68,-0.02))

cx3.legend(frameon=False,loc=(1.1,0.5))
cx3.set_xlim(xmin=0,xmax=1)
cx3.axis('off')
'''
formatter='svg'
#plt.savefig('figures/Fig3_fix.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
