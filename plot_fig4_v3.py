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


#plot position assign
axx=0.1
axy=0.95
bxx=0.1
bxy=0.38
cxx=0.0
cxy=0.0


##vioiline plot
ax=plt.axes((axx,axy,0.5,0.35))
ax.annotate('a',(-0.25,1.07),xycoords='axes fraction',fontweight='bold')
ax.annotate(r'Median of selected mutant frequency in $\Psi(f^*_{k+1}=f|f_k^*)$',(-0.20,1.07),xycoords='axes fraction')
#bx,_=violinplot_from_histogram(bx,ext_pdf[mask],fbins_s[mask],positions=f0[mask],side='left',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C0',alpha=0.5,marker=5,showmeans=False)
#bx,_=violinplot_from_histogram(ax,ext_pdf[mask],fbins_t[mask],positions=f0[mask],side='right',width=f0[mask][1]-f0[mask][0],yoffset=f0[mask],color='C1',alpha=0.5,marker=4,showmeans=False)
#ax.scatter(f0[mask],(sel_median-f0)[mask],c='C0',marker=5)
#ax.scatter(f0[mask],(ext_median-f0)[mask],c='C1',marker=4)
ax.hlines(0,0,1,colors='black',ls=':')
ax.set_xlim(xmin=0,xmax=1)
#ax.set_ylim(ymin=-0.0401,ymax=0.0601)
ax.set_ylim(ymin=-0.022,ymax=0.022)
ax.set_ylabel(r'$f-f^*_k$')
ax.set_xlabel(r" Selected mutant frequency in cycle $k$, $f^*_k$")

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
bx.annotate(r'Distribution of frequency changes in each steps',(-0.20,1.07),xycoords='axes fraction')
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
bx.set_xlabel(r" Selected mutant frequency in cycle $k$, $f^*_k$")

#legend for boxplot
#lcx=plt.axes((0.48,0.1,0.1,0.13,))
#lcx.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
#lcx.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)

#Two strain case
#large frequency
cx1=plt.axes((cxx,cxy,0.4,0.25))
cx1.annotate('c',(-0.,0.8),xycoords='axes fraction',fontweight='bold')
cx1.annotate(r'Success or Fail according to the target composition',(0.05,0.8),xycoords='axes fraction')

#frame
#cx1.hlines(0,0,1,colors='black',lw=3)
cx1.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
#ax.fill_between([0.3,0.7],0.02,0,color='lightgray')
cx1.fill_between([0,1],0.02,0,color='green',alpha=0.3)
#cx1.fill_between([0.3,1],0.02,0,color='red',alpha=0.5)
#ax.annotate('Fail',xy=(0.45,0.005))
cx1.annotate(r'$0$',xy=(0,-0.02))
cx1.annotate(r'$1$',xy=(0.95,-0.02))
#cx1.annotate('Mutant\nfrequency',xy=(1.05,-0.005))
#cx1.annotate('a',xy=(-0.0,0.95),weight='bold',xycoords='axes fraction')
#cx1.annotate(' Two-strain',xy=(0.05,0.95),xycoords='axes fraction')

cx1.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.9
ini1=0.2
h1=0.01
ini2=0.5
h2=0.01
ini3=0.8
h3=0.01
#cx1.scatter(tgt,0.015,c='black',marker='s') #target
cx1.vlines(tgt,0,0.02,colors='black',ls=':',label='Target') #target
#cx1.annotate('Target',xy=(tgt,0.015),xytext=(tgt-0.1,0.035),arrowprops=dict(arrowstyle='->'))
cx1.scatter(ini1,h1,c='black',marker='o',label='Initial') #initial 1
#cx1.annotate('Initial',xy=(ini1,0.015),xytext=(ini1-0.05,0.03),arrowprops=dict(arrowstyle='->'))
cx1.scatter(ini2,h2,c='black',marker='o') #initial 2
#cx1.annotate('Initial',xy=(ini2,0.015),xytext=(ini2-0.05,0.03),arrowprops=dict(arrowstyle='->'))
cx1.scatter(ini3,h3,c='black',marker='o') #initial 2
#cx1.annotate('Initial',xy=(ini3,0.015),xytext=(ini3-0.05,0.03),arrowprops=dict(arrowstyle='->'))
#cx1.annotate('',xy=(ini3,0.01),xytext=(ini3-0.05,0.01),arrowprops=dict(arrowstyle='->'))

cx1.annotate('',xy=(ini2,h2),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
cx1.annotate('',xy=(ini3,h3),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
cx1.annotate('',xy=(tgt,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))
#ax.plot([0.25,0.9],[0.01,0.01],ls=':',color='black')
#ax.scatter([0.5],[0.01],marker='x',color='red')
#cx1.annotate('Fail',xy=(0.5,0.01),xytext=(0.45,-0.03),arrowprops=dict(arrowstyle='->'))
#ax.scatter([0.15],[0.01],facecolor='none',edgecolor='blue')
#cx1.annotate('Success',xy=(0.15,0.01),xytext=(0.05,-0.03),arrowprops=dict(arrowstyle='->'))


cx1.vlines([0.3,0.7],0.005,-0.005,colors='black')
cx1.annotate(r'$f^L$',xy=(0.28,-0.02))
cx1.annotate(r'$f^U$',xy=(0.68,-0.02))

cx1.set_xlim(xmin=0,xmax=1)

#cx1.legend(frameon=False,loc=(0,-0.1))
cx1.axis('off')

#Two strain case
#mid frequency
cx2=plt.axes((cxx,cxy-0.15,0.4,0.25))

#frame
#cx1.hlines(0,0,1,colors='black',lw=3)
cx2.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
#ax.fill_between([0.3,0.7],0.02,0,color='lightgray')
cx2.fill_between([0,1],0.02,0,color='red',alpha=0.3)
#cx1.fill_between([0.3,1],0.02,0,color='red',alpha=0.5)
#ax.annotate('Fail',xy=(0.45,0.005))
cx2.annotate(r'$0$',xy=(0,-0.02))
cx2.annotate(r'$1$',xy=(0.95,-0.02))
#cx2.annotate('Mutant\nfrequency',xy=(1.05,-0.005))
#cx2.annotate('(a) Two-strain',xy=(-0.2,0.95),xycoords='axes fraction')

cx2.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.5
ini1=0.2
h1=0.010
ini2=0.5
h2=0.01
ini3=0.8
h3=0.01
#cx2.scatter(tgt,0.015,c='black',marker='s') #target
cx2.vlines(tgt,0,0.02,colors='black',ls=':',label='Target') #target
#cx2.annotate('Target',xy=(tgt,0.015),xytext=(tgt-0.1,0.035),arrowprops=dict(arrowstyle='->'))
cx2.scatter(ini1,h1,c='black',marker='o',label='Initial(low)') #initial 1
#cx1.annotate('Initial',xy=(ini1,0.015),xytext=(ini1-0.05,0.03),arrowprops=dict(arrowstyle='->'))
cx2.scatter(ini2,h2,c='black',marker='o',label='Initial(mid)') #initial 2
#cx1.annotate('Initial',xy=(ini2,0.015),xytext=(ini2-0.05,0.03),arrowprops=dict(arrowstyle='->'))
cx2.scatter(ini3,h3,c='black',marker='o',label='Initial(high)') #initial 2
#cx1.annotate('Initial',xy=(ini3,0.015),xytext=(ini3-0.05,0.03),arrowprops=dict(arrowstyle='->'))
#cx1.annotate('',xy=(ini3,0.01),xytext=(ini3-0.05,0.01),arrowprops=dict(arrowstyle='->'))

cx2.annotate('',xy=(ini2,h2),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
cx2.annotate('',xy=(0.7,h2),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
cx2.annotate('',xy=(0.7,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))
#ax.plot([0.25,0.9],[0.01,0.01],ls=':',color='black')
#ax.scatter([0.5],[0.01],marker='x',color='red')
#cx1.annotate('Fail',xy=(0.5,0.01),xytext=(0.45,-0.03),arrowprops=dict(arrowstyle='->'))
#ax.scatter([0.15],[0.01],facecolor='none',edgecolor='blue')
#cx1.annotate('Success',xy=(0.15,0.01),xytext=(0.05,-0.03),arrowprops=dict(arrowstyle='->'))
#cx2.annotate('',xy=(0.7,h2),xytext=(tgt,h2),arrowprops=dict(arrowstyle='fancy'))

cx2.vlines([0.3,0.7],0.005,-0.005,colors='black')
cx2.annotate(r'$f^L$',xy=(0.28,-0.02))
cx2.annotate(r'$f^U$',xy=(0.68,-0.02))

cx2.set_xlim(xmin=0,xmax=1)
cx2.axis('off')

#Two strain case
#small frequency
cx3=plt.axes((cxx,cxy-0.3,0.4,0.25))

#frame
#cx1.hlines(0,0,1,colors='black',lw=3)
cx3.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
#ax.fill_between([0.3,0.7],0.02,0,color='lightgray')
cx3.fill_between([0,0.3],0.025,0,color='green',alpha=0.3, label='success')
cx3.fill_between([0.3,1],0.025,0,color='red',alpha=0.3,label='fail')
#ax.annotate('Fail',xy=(0.45,0.005))
cx3.annotate(r'$0$',xy=(0,-0.02))
cx3.annotate(r'$1$',xy=(0.95,-0.02))
cx3.annotate('Mutant frequency',xy=(0.5,0.2),xycoords='axes fraction',ha='center',va='top')
#cx2.annotate('(a) Two-strain',xy=(-0.2,0.95),xycoords='axes fraction')

cx3.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.05
ini1=0.2
h1=0.015
ini2=0.5
h2=0.015
ini3=0.8
h3=0.015
#cx3.scatter(tgt,0.015,c='black',marker='s') #target

#cx3.annotate('Target',xy=(tgt,0.015),xytext=(tgt-0.1,0.035),arrowprops=dict(arrowstyle='->'))
cx3.scatter(ini1,h1,c='black',marker='o',label='Initial') #initial 1
#cx1.annotate('Initial',xy=(ini1,0.015),xytext=(ini1-0.05,0.03),arrowprops=dict(arrowstyle='->'))
cx3.scatter(ini2,h2,c='black',marker='o') #initial 2
#cx1.annotate('Initial',xy=(ini2,0.015),xytext=(ini2-0.05,0.03),arrowprops=dict(arrowstyle='->'))
cx3.scatter(ini3,h3,c='black',marker='o') #initial 2
cx3.vlines(tgt,-0.005,0.025,colors='black',ls=':',label='Target') #target
#cx1.annotate('Initial',xy=(ini3,0.015),xytext=(ini3-0.05,0.03),arrowprops=dict(arrowstyle='->'))
#cx1.annotate('',xy=(ini3,0.01),xytext=(ini3-0.05,0.01),arrowprops=dict(arrowstyle='->'))

cx3.annotate('',xy=(tgt,h1),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
cx3.annotate('',xy=(0.7,h2),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
cx3.annotate('',xy=(0.7,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))
#ax.plot([0.25,0.9],[0.01,0.01],ls=':',color='black')
#ax.scatter([0.5],[0.01],marker='x',color='red')
#cx1.annotate('Fail',xy=(0.5,0.01),xytext=(0.45,-0.03),arrowprops=dict(arrowstyle='->'))
#ax.scatter([0.15],[0.01],facecolor='none',edgecolor='blue')
#cx1.annotate('Success',xy=(0.15,0.01),xytext=(0.05,-0.03),arrowprops=dict(arrowstyle='->'))
#cx3.annotate('',xy=(0.7,h2),xytext=(tgt,h2),arrowprops=dict(arrowstyle='fancy'))

cx3.vlines([0.3,0.7],0.005,-0.005,colors='black')
cx3.annotate(r'$f^L$',xy=(0.28,-0.02))
cx3.annotate(r'$f^U$',xy=(0.68,-0.02))

cx3.legend(frameon=False,loc=(1.1,0.5))
cx3.set_xlim(xmin=0,xmax=1)
cx3.axis('off')

formatter='svg'
#plt.savefig('Fig4_v3.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
#plt.savefig('Fig4_v2.png',dpi=300,bbox_inches='tight')
plt.show()
