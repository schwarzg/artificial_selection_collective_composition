import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "arial"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.rcParams["mathtext.fontset"]='stix'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'STIXGeneral:italic:bold'
plt.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
plt.rcParams['font.size']=12

from analytic_results import *

mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
rhat=0.05
nens=300
ncycle=10
tcycle=np.log(ncomm+1)/r

def get_f(w,m):
    return np.divide(m,w+m)

barN02=[]
barf=[]
sigw2=[]
sigm2=[]
sigf2=[]
N0s= np.arange(500,5000,100)[::-1]
for N0 in N0s:
    m0=int(N0*0.5)#np.clip(np.random.poisson(200,size=nens),0,1000)
    w0=N0-m0
    c0=get_f(w0,m0)
    barw=barw_th(np.mean(w0),tcycle,r,mu)    
    sig2w=sig2w_th(np.mean(w0),np.var(w0),tcycle,r,mu)
    barm=barm_th(np.mean(m0),np.mean(w0),tcycle,r,mu,s)
    sig2m=sig2m_th(np.mean(m0),np.var(m0),np.mean(w0),np.var(w0),-np.var(m0),tcycle,r,mu,s)    
    sig2c=sig2c_th_v2(np.mean(c0),np.var(c0),N0,tcycle,r,mu,s)

    barN02.append((barw+barm)**2)
    barf.append(barm/(barw+barm))
    sigw2.append(sig2w)
    sigm2.append(sig2m)
    sigf2.append(sig2c)


fig,ax=plt.subplots(1,3,figsize=(12,3))


#ax=plt.axes((0.1,0.1,0.3,0.3))
#ax[0]=plt.axes((0.53,0.1,0.3,0.3))
'''
ax.plot(N0s,sigw2,label=r'$\sigma_w^2(\tau)$')    
ax.plot(N0s,sigm2,label=r'$\sigma_m^2(\tau)$')
ax.plot([1500,3000],[45000,90000],ls=':',c='black')
ax.annotate(r'$\sim N_0$',(2000,45000))
ax.set_xlabel(r'Newborn collective size $N_0$')    
ax.set_ylabel(r'Variances')    
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(frameon=False,handlelength=0.7)
ax.set_xticks([500,1000,5000])
ax.set_xticklabels([500,1000,5000])
ax[0].plot(N0s,barN02)    
ax[0].plot([1500,3000],[150000000,600000000],ls=':',c='black')
ax[0].annotate(r'$\sim N_0^2$',(2000,150000000))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel(r'Newborn collective size $N_0$')    
ax[0].set_ylabel(r'$\overline{N}(\tau)^2$')    
ax[0].set_xticks([500,1000,5000])
ax[0].set_xticklabels([500,1000,5000])
'''
#ax[0]=plt.axes((0.1,0.1,0.3,0.3))
ax[0].annotate('a',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')

from analytic_results import barc_th_v2,sig2c_th_v2

f0s = np.arange(0,1+0.01*0.5,0.01)
barcs=[]
sig2cs=[]
for f0 in f0s:
    barc=barc_th_v2(f0,N0,tcycle,r,mu,s)
    barcs.append(barc)
    #sig2c=sig2c_th_v2(f0,f0*(1-f0)/N0,N0,tcycle,r,mu,s)
    sig2c=sig2c_th_v2(f0,0,N0,tcycle,r,mu,s)
    sig2cs.append(sig2c)

ax[0].plot(f0s,barcs-f0s,label=r"$\bar{f}_{k+1}-f^*_k$")

#ax[0].set_xlabel(r"Selected F frequency  $f^*_k$")
ax[0].set_xlabel(r"Newborn frequency $f_0$")
#ax[0].set_ylabel(r"$\bar{f}_{k+1}-f^*_k$")
#ax[0].set_ylabel(r"$\sigma_{{f}_{k+1}}$")
ax[0].tick_params(axis='y',labelcolor='C0')
#ax[0].set_ylabel(r"Mean frequency difference",color='C0')
ax[0].set_ylabel(r'$\bar{f}\ \ -f_0$ (Average difference in $f$\nbetween Adults and a fixed Newborn'.replace(r'\n','\n'),color='C0')
#(ymin,ymax)=ax[0].get_ylim()

ax0t=ax[0].twinx()
ax0t.plot(f0s,sig2cs,label=r"$\sigma^2_{{f}_{k+1}}$",c='C1')
ax0t.tick_params(axis='y',labelcolor='C1')
#ax0t.set_ylabel(r"Standard deviation",color='C1')
ax0t.set_ylabel(r'Variance in Adults $f$',color='C1')
#ax0t.set_ylim(ymin=ymin,ymax=ymax)
ax0t.set_ylim(ymax=0.0008)
#ax[0].legend(frameon=False)

#ax[1]=plt.axes((0.68,0.1,0.3,0.3))
ax[1].annotate('b',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')
ax[1].plot(1/N0s,sigf2,c='C1')    
ax[1].plot([1/1000,1/2000],[0.00016,0.00008],ls=':',c='black')
ax[1].annotate(r'$\sim 1/N_0$',(1/1000,0.00012))
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel(r'Newborn size $1/N_0$')    
#ax[1].set_ylabel(r'$\sigma_f^2(\tau)$')    
ax[1].set_ylabel(r'Variance in Adult $f$ ($\sigma_{f,}^2\ \ )$',color='C1')    
ax[1].tick_params(axis='y',labelcolor='C1')
ax[1].set_xticks([1/1000,1/5000])
ax[1].set_xticklabels(['1/1000','1/5000'])
ax[1].set_ylim(ymin=10**-5)
'''
ax[1]t=ax[1].twinx()
ax[1]t.plot(1/N0s,barf)    
ax[1]t.tick_params(axis='y',labelcolor='C1')
#ax[0]t.set_ylabel(r"Standard deviation",color='C1')
ax[1]t.set_ylabel(r'$\sigma^2_{f,\tau}$',color='C1')
#ax[0]t.set_ylim(ymin=ymin,ymax=ymax)
#ax[0]t.set_ylim(ymax=0.0006)
#ax[0].legend(frameon=False)
'''

f0=0.5
N0=1000
ts=np.arange(0.4*tcycle,14*tcycle,0.2*tcycle)
barcs=[]
barws=[]
barms=[]
sig2ms=[]
sig2ws=[]
sig2cs=[]
barcs2=[]
barws2=[]
barms2=[]
sig2ms2=[]
sig2ws2=[]
sig2cs2=[]
for t in ts:
    #formula type meanfield - analytic
    Rt=np.exp(r*t)
    Wt=np.exp(s*t)
    Mt=np.exp(-mu*t)

    num=f0*Wt+mu/s*(1-f0)*(Wt-1)
    den=1-f0+num
    barc=num/den 
    barcs.append(barc)

    num1=f0*Wt*((2-2*f0+2*f0**2)*Rt*Wt-(1-f0)-f0*Wt)
    num2=mu/s*(1-f0)*((1-f0)*((2*s/(r+2*s)-2*f0)*Rt*Wt**2-Wt+r/(r+2*s)+2*f0*Rt*Wt)+2*f0*Wt*(Wt-1)*((1+f0)*Rt-1))

    den=N0*Rt*(1-f0+f0*Wt+mu/s*(1-f0)*(Wt-1))**4

    sig2c=(1-f0)*(num1+num2)/den
    sig2cs.append(sig2c) 
  
    '''
    Rt=np.exp(r*t)
    Wt=np.exp(s*t)
    num=f0*(1-f0)*(2-2*f0+2*f0**2-(1-f0)/Rt/Wt-f0/Rt)
    den=N0*Wt**2*((1-f0)/Wt+f0)**4
    sig2c=num/den
    ''' 
    
    barc=barc_th_v2(f0,N0,t,r,mu,s)
    barcs2.append(barc)
    mu=0
    m0=f0*N0
    w0=(1-f0)*N0
    sig2w0=N0*f0*(1-f0)#(N0**2)*sig2f0/(2*f0**2-2*f0+1)
    sig2m0=sig2w0
    covmw0=-sig2w0
    #print(w0,m0,sig2w0,sig2m0)
    
    w=barw_th(w0,t,r,mu)
    m=barm_th(m0,w0,t,r,mu,s)
    #w= w0*Rt*Mt
    erst=np.exp((r+s)*t)

    N=(w+m)
    N2=N**2
    c=m/N
    sig2w=sig2w_th(w0,sig2w0,t,r,mu)
    sig2m=sig2m_th(m0,sig2m0,w0,sig2m0,-sig2m0,t,r,mu,s)

    sig2c=((1-c**2)*sig2m+c**2*sig2w)/N2
    sig2c=sig2c_th_v2(f0,(1-f0)*f0/N0,N0,t,r,mu,s)
    #sig2c=sig2c_th((1-f0)*N0,f0*(1-f0)*N0,f0*N0,f0*(1-f0)*N0,-f0*(1-f0)*N0,t,r,mu,s)
    #sig2c=sig2c_th_v2(f0,0,N0,t,r,mu,s)
    sig2cs2.append(sig2c) #composite formula type 2

ax[2].annotate('c',(-0.1,1.05),xycoords='axes fraction',fontweight='bold')
#ax[2].plot(ts,np.array(barcs)-f0,c='C1')    
#ax[2].plot(ts,barcs,c='C1')    
#ax[2].plot(ts,barcs2,c='C2')    
print(ts,sig2cs,sig2cs2)
ax[2].plot(ts,sig2cs2,c='C1')    
#ax[2].plot(ts,sig2cs,c='C2',ls=':')    
#ax[2].plot(ts,sig2cs3,c='C3')    
#ax[2].set_xscale('log')
#ax[2].set_yscale('log')
ax[2].set_xlabel(r'Maturation time $\tau$')    
#ax[1].set_ylabel(r'$\sigma_f^2(\tau)$')    
ax[2].set_ylabel(r'Variance in Adult $f$ ($\sigma_{f,}^2\ \ )$',color='C1')    
ax[2].tick_params(axis='y',labelcolor='C1')
#ax[2].tick_params(axis='y',labelcolor='C1')
#ax[2].set_xticks([1/1000,1/5000])
#ax[2].set_xticklabels(['1/1000','1/5000'])
#ax[2].set_ylim(ymin=10**-5)

fig.tight_layout()

formatter='svg'
#plt.savefig('figures/FigS6_v2.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
