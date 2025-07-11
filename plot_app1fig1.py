import numpy as np
from analytic_results import *


import matplotlib.pyplot as plt

#plt.rcParams['text.usetex']=True
plt.rcParams["font.family"] = "arial"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.size"] = 11
plt.rcParams["axes.labelweight"] = "bold"
#plt.rcParams["mathtext.fontset"]='stix'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'STIXGeneral:italic:bold'
plt.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'

from model import model
from tau_leaping import tau_leaping
from growth_alter import *
import scipy.stats as st

mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
rhat=0.05
nens=3000
ncycle=10
tcycle=np.log(3*ncomm+1)/r

proFunc=[
lambda x: r*x[0],
lambda x: mu*x[0],
lambda x: (r+s)*x[1]
]
changeVec=np.array([[1,0],[-1,1],[0,1]]).astype(int)
rxnOrder=np.array([[1,0],[1,0],[0,1]]).astype(int)

growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)

solver=tau_leaping(model=growthmodel,eps=0.01)


N0=1000

fig,ax=plt.subplots(3,3,figsize=(8,6))
ax[0,0].annotate(r'$(S_0,F_0)$',xy=(-0.5,1.1),xycoords='axes fraction')
    
for i,f0 in enumerate([0.01,0.5,0.99]):
    #c0=0.05
    
    m0=np.ones(nens)*int(N0*f0)#np.clip(np.random.poisson(200,size=nens),0,1000)
    w0=N0-m0
    
    barw=barw_th(np.mean(w0),tcycle,r,mu)    
    sig2w=sig2w_th(np.mean(w0),np.var(w0),tcycle,r,mu)
    barm=barm_th(np.mean(m0),np.mean(w0),tcycle,r,mu,s)
    sig2m=sig2m_th(np.mean(m0),np.var(m0),np.mean(w0),np.var(w0),-np.var(m0),tcycle,r,mu,s)    
    def get_f(w,m):
        return np.divide(m,w+m)
    c0=get_f(w0,m0)
    barc=barc_th_v2(np.mean(c0),N0,tcycle,r,mu,s)
    sig2c=sig2c_th_v2(np.mean(c0),np.var(c0),N0,tcycle,r,mu,s)
      
    ws=np.zeros(nens)
    ms=np.zeros(nens)
    ws2=np.zeros(nens)
    ms2=np.zeros(nens)
    from tqdm import tqdm
    for e in tqdm(range(nens)):
        #print(e)
        _,X=solver.run(np.array([w0[e],m0[e]]),0,tcycle)
        ws[e]=X[0,-1]    
        ms[e]=X[1,-1]    
        #ws2[e]=growth_sampling_w(np.mean(w0),np.var(w0),tcycle,r,mu,s)
        #ms2[e]=growth_sampling_m(np.mean(m0),np.var(m0),np.mean(w0),np.var(w0),-np.var(m0),tcycle,r,mu,s)
    cs=get_f(ws,ms)    
    #ax=plt.axes((0.1,0.5,0.4,0.3))
    #ax[0,0]=plt.axes((0.1,0.5,0.30,0.3))
    #ax.hist(w0,bins=10,density=True)
    ax[0,i].annotate(r'$(%d,%d)$'%(w0[0],m0[0]),xy=(0.5,1.1),xycoords='axes fraction',ha='center')
    ax[0,i].hist(ws,bins=30,density=True,alpha=0.4,label='tau')
    #ax[0,i].hist(ws2,bins=30,density=True,alpha=0.4,label='samp')
    width=max(abs(barw-max(ws)),abs(barw-min(ws)))
    wrange=np.linspace(max(0,barw-width),barw+width,100)    
    pdfw=st.norm.pdf(wrange,loc=barw,scale=np.sqrt(sig2w))
    ax[0,i].plot(wrange,pdfw,label='Gauss')
    ax[0,i].set_xlabel(r'$S$')
    ax[0,i].set_yscale('log')
    ax[0,i].legend(frameon=False,handlelength=0.5)
    
    #bx=plt.axes((0.1,0.1,0.4,0.3))
    #ax[1,0]=plt.axes((0.55,0.5,0.30,0.3))
    #bx.hist(m0,bins=10,density=True)
    ax[1,i].hist(ms,bins=30,density=True,alpha=0.4)        
    #ax[1,i].hist(ms2,bins=30,density=True,alpha=0.4)        
    width=max(abs(barm-max(ms)),abs(barm-min(ms)))
    mrange=np.linspace(max(0,barm-width),barm+width,100)    
    pdfm=st.norm.pdf(mrange,loc=barm,scale=np.sqrt(sig2m))
    ax[1,i].plot(mrange,pdfm)
    ax[1,i].set_xlabel(r'$F$')
    ax[1,i].set_yscale('log')
    
    #ax[2,i]=plt.axes((0.1,0.1,0.35,0.3))
    ax[2,i].hist(cs,density=True,alpha=0.4,bins=30)
    #print(barc,sig2c,'?=',np.mean(cs),np.var(cs))    
    def func(f0,t):
        Rt=np.exp(r*t)
        Wt=np.exp(s*t)
        num1=(1-f0)**2 * (f0*(1-f0)*Rt*Wt**2+f0*Wt*(Rt*Wt*-1)+mu/s*(1-f0)*((2*s/(r+2*s)-2*f0)*Rt*Wt**2+2*f0*Rt*Wt-Wt+r/(r+2*s)))
        num2=(f0*Wt+mu/s*(1-f0)*(Wt-1))**2 * (f0*(1-f0)*Rt+(1-f0)*(Rt-1))
        den=N0*Rt*(1-f0+f0*Wt+mu/s*(1-f0)*(Wt-1))**4
        return (num1+num2)/den
    #sig2c=func(np.mean(c0),tcycle)
    width=max(abs(barc-max(cs)),abs(barc-min(cs)))
    crange=np.linspace(max(0,barc-width),min(1,barc+width),100)    
    pdfc=st.norm.pdf(crange,loc=barc,scale=np.sqrt(sig2c))
    ax[2,i].plot(crange,pdfc)
    ax[2,i].set_xlabel(r'$f$')
    ax[2,i].set_yscale('log')

ax[0,0].set_ylabel(r'$P(S)$')
ax[1,0].set_ylabel(r'$P(F)$')
ax[2,0].set_ylabel(r'$P(f)$')
    
fig.tight_layout()
#plt.savefig('figures/FigS1_v2.svg',dpi=300,bbox_inches='tight',format='svg')
plt.show()
