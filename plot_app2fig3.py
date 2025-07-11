'''

Test stability of median and mean in farely skewed distribution 

'''


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'DejaVu Sans:italic:bold'
plt.rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'

import time

import scipy.stats as st
import scipy.integrate as intg
from scipy.integrate import simpson
import scipy.special as spc
from analytic_results import *

from one_step_functions import *
from model import model
from tau_leaping import tau_leaping
from direct import direct
from growth_alter import *
from selection import select
from selection import community_function as cf
from selection import score_function as sf
from reproduction import hypergeometric_reproduce
from stats_custom_tools import * 

import matplotlib.pyplot as plt
from tqdm import tqdm

#parameter prepare
mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=1000
rhat=0 if s>0 else 1
nens=1000
ncycle=10

ntest=1#00


f0=0.5

fig,ax=plt.subplots(1,3,figsize=(12,3))

labels=['a','b','c']
ncomms=[10,100,1000]    
for i in range(3):
    label=labels[i]
    ncomm=ncomms[i]
    tcycle=np.log(1000)/r
    ############################
    # Evaluation part
    ############################    
    fbins=np.linspace(np.maximum(0,f0-0.25),np.minimum(1,f0+0.25),100)
    #fbins=np.linspace(np.maximum(0,f0-0.05),np.minimum(1,f0+0.05),30)
    skewu=0.1 if s>0 else 0
    skewd=0.1 if s<0 else 0
    foribins=np.linspace(np.maximum(0,f0-0.15-skewd),np.minimum(1,f0+0.15+skewu),100)     #extend the range to describe the distribution more
    pdfforplot=transition_probability_iid_pdf(foribins[:-1],f0,tcycle,N0,r,mu,s)
    #plt.plot(foribins[:-1],pdfforplot)
    #plt.hist(dats,bins=foribins,density=True)
    pdf=transition_probability_iid_pdf(fbins[:-1],f0,tcycle,N0,r,mu,s)
    pdfnorm=simpson(pdf,dx=fbins[1]-fbins[0])
    if f0-0.15<0 or 1<f0+0.15:
        pdf=pdf/pdfnorm        #if pdf range is not enough, reweight the pdf for sure completeness
    cdf=np.cumsum(pdf*(fbins[1]-fbins[0]))
    print(cdf[-1])
    extpdf=ncomm*(1-cdf)**(ncomm-1)*pdf #Same dimension with dat
    extpdf=extpdf/simpson(extpdf,dx=(fbins[1]-fbins[0]))
    median=quantile(0.5,fbins,extpdf,x0=np.min(fbins))
    mean=simpson(fbins[:-1]*extpdf,dx=(fbins[1]-fbins[0]))
    print(median-mean)
    #plt.plot(fbins[:-1],extpdf)
    #pdf=transition_probability_iid_pdf(fbins[:-1],f0,tcycle,N0,r,mu,s)
    #pdfnorm=simpson(pdf,dx=fbins[1]-fbins[0])
    #if f0-0.15<0 or 1<f0+0.15:
    #    pdf=pdf/pdfnorm        #if pdf range is not enough, reweight the pdf for sure completeness
    #plt.show()
    
    ############################
    # Simulation part
    ############################
    #growth phase
    #prepare data space
    w_sel=np.zeros(nens)
    m_sel=np.zeros(nens)
    fsampmean=np.array([])
    fsampstd=np.array([])
    #ensmble loop
    dats=[]
    for e in tqdm(range(nens),leave=False):
        #print("ensemble ",e)
        #reproduction phase    
        m0sel=np.random.binomial(N0,f0,size=ncomm)    
        m0sel=np.sort(m0sel)
        w0sel=N0-m0sel
        
        #For each community
        lastw=np.zeros(ncomm)
        lastm=np.zeros(ncomm)
        for j in range(ncomm):
            #Growth phase
            wsamp=growth_sampling_w(w0sel[j],0,tcycle,r,mu,s)
            msamp=growth_sampling_m(m0sel[j],0,w0sel[j],0,0,tcycle,r,mu,s)
    
            #store for selection
            lastw[j]=wsamp
            try:
                lastm[j]=msamp
            except:
                print(msamp)
            
            #get and store frequency
            dats.append(cf(wsamp,msamp))
        
        #Selection phase
        ind_sel=select(lastw,lastm,r=r,s=s,rhat=rhat)
        w_sel[e]=lastw[ind_sel]
        m_sel[e]=lastm[ind_sel]
        
    cfs_sel=cf(w_sel,m_sel)
    dmedians=cfs_sel-median
    dmeans=cfs_sel-mean
    
    #ax[i].hist(dmedians,bins=30,label='median',alpha=0.5)
    #ax[i].hist(dmeans,bins=30,label='mean',alpha=0.5)
    p,_,_=ax[i].hist(cfs_sel,bins=30,alpha=0.5)
    ax[i].vlines(np.median(cfs_sel),0,max(p),colors='C1')
    #ax[i].vlines(median,0,max(p),colors='C1')
    #ax[i].vlines(mean,0,max(p),colors='C4')
    ax[i].vlines(np.mean(cfs_sel),0,max(p),colors='C4')
    ax[i].annotate(label,xy=(-0.1,1.1),xycoords='axes fraction')
ax[0].set_xlabel(r'F frequency $f$')
ax[1].set_xlabel(r'F frequency $f$')
ax[2].set_xlabel(r'F frequency $f$')
ax[0].set_ylabel('Probability density')
ax[0].legend(frameon=False)
fig.tight_layout()
#plt.savefig('figures/figS6.52.svg',format='svg',dpi=300,bbox_inches='tight')
plt.show()
