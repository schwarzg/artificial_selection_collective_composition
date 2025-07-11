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
Nf=(ncomm+2)*N0        #Nf is set to be (ncomm+1)*N0
fbins=np.linspace(0,1,100)

#get data

#frequency loop
fls=np.array([])
fle=np.array([])
fus=np.array([])
fue=np.array([])
nens=300
for ncomm in [4,6,8,10,20,40,60,80,100]:
    folder="data/cond/conditional_probability_fixtau_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
    foldere="data/cond/conditional_probability_fixtau_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix"%(N0,r,mu,s,ncomm,1000)
    seldat=np.loadtxt(folder+"_smedian.dat")
    extdat=np.loadtxt(foldere+"_emedian.dat")[1:]
    if ncomm==10 or ncomm==20 or ncomm==40:
        seldat=seldat[2:]#np.loadtxt(folder+"_smedian.dat")
        #extdat=extdat[1:]#np.loadtxt(folder+"_emedian.dat")
    f0=np.arange(0.02,0.99,0.01)
    
    hatsel=savgol_filter(seldat-f0,11,2)    #smoothing
    selmean_interp=itp.interp1d(f0,hatsel)
    extmean_interp=itp.interp1d(f0,extdat-f0)
    fl_sel=None
    lstart=0.0
    while fl_sel is None:
        try:
            fl_sel=opt.root(selmean_interp,lstart).x
        except:
            lstart=lstart+0.05
            continue
        
    lstart=0.99
    fu_sel=None
    while fu_sel is None:
        try:
            fu_sel=opt.root(selmean_interp,lstart).x
        except:
            lstart=lstart-0.05
            continue
    fl_ext=None
    lstart=0.0
    while fl_ext is None:
        try:
            fl_ext=opt.root(extmean_interp,lstart).x
        except:
            lstart=lstart+0.05
            continue
    lstart=0.99
    fu_ext=None
    while fu_ext is None:
        try:
            fu_ext=opt.root(extmean_interp,lstart).x
        except:
            lstart=lstart-0.05
            continue
        if fu_ext < 0.5:
            fu_ext=1
    print(ncomm,fl_sel,fu_sel,fl_ext,fu_ext)
    fls=np.append(fls,fl_sel)
    fus=np.append(fus,fu_sel)
    fle=np.append(fle,fl_ext)
    fue=np.append(fue,fu_ext)
#N0 data
nens=1000
ncomm=10
fls2=np.array([])
fle2=np.array([])
fus2=np.array([])
fue2=np.array([])
N0s=[700,800,900,1000,2000,4000,6000]
for N0 in N0s:
    folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
    foldere="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix"%(N0,r,mu,s,ncomm,nens)
    seldat=np.loadtxt(folder+"_smedian.dat")[1:]
    extdat=np.loadtxt(foldere+"_emedian.dat")
    f0=np.arange(0.01,0.99,0.01)#[:len(seldat)]
    
    hatsel=savgol_filter(seldat-f0,11,2)    #smoothing
    selmean_interp=itp.interp1d(f0,hatsel,fill_value='extrapolate')
    extmean_interp=itp.interp1d(f0,extdat-f0,fill_value='extrapolate')
    
    #search range when fl is approachable
    left1=0.01
    right1=0.02
    while selmean_interp(left1)*selmean_interp(right1)>0 and right1<1.0:     #until the sign is different
        right1=right1+0.01
    left2=0.51
    right2=0.52
    while selmean_interp(left2)*selmean_interp(right2)>0 and right2<1.0:     #until the sign is different
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
    while extmean_interp(left1)*extmean_interp(right1)>0 and right1<1.0:     #until the sign is different
        right1=right1+0.01
    left2=0.51
    right2=0.52
    while extmean_interp(left2)*extmean_interp(right2)>0 and right2<1.0:     #until the sign is different
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

ncomm=10
N0=1000
fls3=np.array([])
fle3=np.array([])
fus3=np.array([])
fue3=np.array([])
logts=np.append(range(11,21,1),range(21,101,10))#[,800,900,1000,2000,4000,6000]
for logt in logts:
    folder="data/cond/conditional_probability_taulog%s_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(logt,N0,r,mu,s,ncomm,nens)
    foldere=folder#"data/cond/conditional_probability_taulog%s_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix"%(logt,ncomm,r,mu,s,ncomm,nens)
    seldat=np.loadtxt(folder+"_smedian.dat")
    extdat=np.loadtxt(foldere+"_emedian.dat")
    f0=np.arange(0.01,0.99,0.01)#[:len(seldat)]
    
    hatsel=savgol_filter(seldat-f0,11,2)    #smoothing
    selmean_interp=itp.interp1d(f0,hatsel,fill_value='extrapolate')
    extmean_interp=itp.interp1d(f0,extdat-f0,fill_value='extrapolate')
    plt.plot(f0,extdat-f0)
    #plt.show()
    
    #search range when fl is approachable
    left1=0.01
    right1=0.02
    while selmean_interp(left1)*selmean_interp(right1)>0 and right1<1.0:     #until the sign is different
        right1=right1+0.01
    left2=0.51
    right2=0.52
    while selmean_interp(left2)*selmean_interp(right2)>0 and right2<1.0:     #until the sign is different
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
    while extmean_interp(left1)*extmean_interp(right1)>0 and right1<1.0:     #until the sign is different
        right1=right1+0.01
    left2=0.51
    right2=0.52
    while extmean_interp(left2)*extmean_interp(right2)>0 and right2<1.0:     #until the sign is different
        right2=right2+0.01
    print(right1,right2)
    if np.abs(right1-1)<0.1:
        fl_ext=0.5
        fu_ext=0.5
    else:
        fl_ext=opt.root(extmean_interp,right1+0.1).x
        fu_ext=opt.root(extmean_interp,right2).x
    print(logt,fl_sel,fu_sel,fl_ext,fu_ext)
    fls3=np.append(fls3,fl_sel)
    fus3=np.append(fus3,fu_sel)
    fle3=np.append(fle3,fl_ext)
    fue3=np.append(fue3,fu_ext)
plt.show()
'''
#N0 data
nens=1000
ncomm=10
N0=1000
fls3=np.array([])
fle3=np.array([])
fus3=np.array([])
fue3=np.array([])
ts=np.arange(0,10.1,0.5)
for dt in ts:
    if dt==0:
        folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
    else:
        folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_tau%.2f_nens%s"%(N0,r,mu,s,ncomm,tcycle+dt,nens)
    seldat=np.loadtxt(folder+"_smedian.dat")
    extdat=np.loadtxt(folder+"_emedian.dat")
    f0=np.arange(0.00,0.99,0.01)[:len(seldat)]
    
    hatsel=savgol_filter(seldat-f0,11,2)    #smoothing
    selmean_interp=itp.interp1d(f0,hatsel,fill_value='extrapolate')
    extmean_interp=itp.interp1d(f0,extdat-f0,fill_value='extrapolate')
    
    #search range when fl is approachable
    left1=0.01
    right1=0.02
    while selmean_interp(left1)*selmean_interp(right1)>0 and right1<1.0:     #until the sign is different
        right1=right1+0.01
    left2=0.51
    right2=0.52
    while selmean_interp(left2)*selmean_interp(right2)>0 and right2<1.0:     #until the sign is different
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
    while extmean_interp(left1)*extmean_interp(right1)>0 and right1<1.0:     #until the sign is different
        right1=right1+0.01
    left2=0.51
    right2=0.52
    while extmean_interp(left2)*extmean_interp(right2)>0 and right2<1.0:     #until the sign is different
        right2=right2+0.01
    print(right1,right2)
    if np.abs(right1-1)<0.1:
        fl_ext=0.5
        fu_ext=0.5
    else:
        fl_ext=opt.root(extmean_interp,right1).x
        fu_ext=opt.root(extmean_interp,right2).x
    print(N0,fl_sel,fu_sel,fl_ext,fu_ext)
    fls3=np.append(fls2,fl_sel)
    fus3=np.append(fus2,fu_sel)
    fle3=np.append(fle2,fl_ext)
    fue3=np.append(fue2,fu_ext)
'''
#draw
shadeup=1.1*np.ones(len(N0s))
shadedw=np.zeros(len(N0s))
print(shadeup,shadedw,fue2)
ax=plt.axes((0.1,0.1,0.3,0.4))
ax.annotate('a',(-0.22,1.08),xycoords='axes fraction',fontweight='bold')
ax.annotate('Success',(0.05,0.85),xycoords='axes fraction')
ax.annotate('Fail',(0.6,0.75),xycoords='axes fraction')
ax.set_ylim(ymin=0,ymax=1)
ax.plot(N0s,fus2,c='C0',label='Sim,$f^U$',marker='v',ms=4)
ax.plot(N0s,fue2,c='C1',label='Th,$f^U$',marker='^',ms=4)
ax.fill_between(N0s,shadeup,fue2,color='lightgray')
ax.plot(N0s,fls2,c='C0',label='Sim,$f^L$',marker='v',ms=4,ls='--')
ax.plot(N0s,fle2,c='C1',label='Th,$f^L$',marker='^',ms=4,ls='--')
ax.fill_between(N0s,fle2,shadedw,color='lightgray')
ax.set_xscale('log')
ax.set_xlim(xmin=300,xmax=6000)
ax.fill_between([1,700],[1.1,1.1],color='lightgray')
ax.set_ylabel(r'Target Frequency $\hat{f}$')
ax.set_xlabel(r'Newborn collective size $N_0$')
ax.set_xticks([500,1000,5000])
ax.set_xticklabels([500,1000,5000])
ax.legend(frameon=False)

def N0_allsuccess():
    b=0.5
    R=np.exp(r*tcycle)
    W=np.exp(s*tcycle)
    E=-1.48
    bunja=E**2*(b*R*W**2+W*(R*W-1)+mu/s*(2*s/(r+2*s)*R*W**2-W+r/(r+2*s))+(W+mu/s*(W-1)**2)*(b*R+R-1))
    bunmo=R*b**3*(1+W+mu/s*(W-1))**2*(1-(W+mu/s*(W-1)))**2
    return bunja/bunmo
def N0_allsuccess2():
    R=np.exp(r*tcycle)
    W=np.exp(s*tcycle)
    E=-1.48
    bunja=8*(W*(3*R*W-W-1)+mu/s*W*(((2*s/(r+2*s)-1)*R*W-1+R+r/(r+2*s))+(3*R-2)*(W-1))) #E**2*(b*R*W**2+W*(R*W-1)+mu/s*(2*s/(r+2*s)*R*W**2-W+r/(r+2*s))+(W+mu/s*(W-1)**2)*(b*R+R-1))
    bunmo=R*(1-(W+mu/s*(W-1))**2)**2
    return E**2*bunja/bunmo

N0bound=N0_allsuccess()
print(N0bound)
N0bound=N0_allsuccess2()
print(N0bound)

ax.vlines(N0bound,0,1,colors='C2',ls=':')
ax.annotate(r'$\breve{N_0}$',(0.3,0.5),xycoords='axes fraction')


bx=plt.axes((0.5,0.1,0.3,0.4))
fs=[4,6,8,10,20,40,60,80,100]
shadeup=1.1*np.ones(len(fs))
shadedw=np.zeros(len(fs))
bx.annotate('b',(-0.22,1.08),xycoords='axes fraction',fontweight='bold')
bx.annotate('Success',(0.65,0.93),xycoords='axes fraction')
bx.annotate('Fail',(0.05,0.45),xycoords='axes fraction')
bx.set_ylim(ymin=0,ymax=1)
bx.set_xscale('log')
bx.plot(fs,fus,c='C0',label='Sim,$f^U$',marker='v',ms=4)
bx.plot(fs,fue,c='C1',label='Th,$f^U$',marker='^',ms=4)
bx.fill_between(fs,shadeup,fue,color='lightgray')
bx.plot(fs,fls,c='C0',label='Sim,$f^L$',marker='v',ms=4,ls='--')
bx.plot(fs,fle,c='C1',label='Th,$f^L$',marker='^',ms=4,ls='--')
bx.fill_between(fs,fle,shadedw,color='lightgray')
bx.set_xlabel(r'Number of collectives $g$')
bx.set_xticks([5,10,50,100])
bx.set_xticklabels([5,10,50,100])
bx.set_xlim(4,100)
'''
cx=plt.axes((1.1,0.1,0.3,0.4))
fs=dt#[4,6,8,10,20,40,60,80,100]
shadeup=1.1*np.ones(len(fs))
shadedw=np.zeros(len(fs))
cx.annotate('b',(-0.22,1.08),xycoords='axes fraction',fontweight='bold')
cx.annotate('Success',(0.65,0.93),xycoords='axes fraction')
cx.annotate('Fail',(0.05,0.45),xycoords='axes fraction')
cx.set_ylim(ymin=0,ymax=1)
cx.set_xscale('log')
cx.plot(fs,fus3,c='C0',label='Sim,$f^U$',marker='v',ms=4)
cx.plot(fs,fue3,c='C1',label='Th,$f^U$',marker='^',ms=4)
cx.fill_between(fs,shadeup,fue,color='lightgray')
cx.plot(fs,fls3,c='C0',label='Sim,$f^L$',marker='v',ms=4,ls='--')
cx.plot(fs,fle3,c='C1',label='Th,$f^L$',marker='^',ms=4,ls='--')
cx.fill_between(fs,fle,shadedw,color='lightgray')
cx.set_xlabel(r'Number of collectives $g$')
cx.set_xticks([5,10,50,100])
cx.set_xticklabels([5,10,50,100])
cx.set_xlim(4,100)

'''
cx=plt.axes((0.9,0.1,0.3,0.4))
#logts=[4,6,8,10,20,40,60,80,100]
ts=np.log(np.array(logts))/r
shadeup=1.1*np.ones(len(logts))
shadedw=np.zeros(len(logts))
cx.annotate('c',(-0.22,1.08),xycoords='axes fraction',fontweight='bold')
cx.annotate('Success',(0.65,0.93),xycoords='axes fraction')
cx.annotate('Fail',(0.05,0.45),xycoords='axes fraction')
cx.plot(ts,fus3,c='C0',label='Sim,$f^U$',marker='v',ms=4)
cx.plot(ts,fue3,c='C1',label='Th,$f^U$',marker='^',ms=4)
cx.fill_between(ts,shadeup,fue3,color='lightgray')
cx.plot(ts,fls3,c='C0',label='Sim,$f^L$',marker='v',ms=4,ls='--')
cx.plot(ts,fle3,c='C1',label='Th,$f^L$',marker='^',ms=4,ls='--')
cx.fill_between(ts,fle3,shadedw,color='lightgray')
cx.fill_between(ts,fue3,fle3,color='gold')
cx.set_xlabel(r'Maturation time $\tau$')
cx.set_ylim(ymin=0,ymax=1)
cx.set_xlim(ts[0],ts[-1])
#cx.set_xscale('log')
#cx.set_xticks([5,10,50,100])
#cx.set_xticklabels([5,10,50,100])
formatter='svg' #or png
#plt.savefig('figures/Fig4_fix2.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()

