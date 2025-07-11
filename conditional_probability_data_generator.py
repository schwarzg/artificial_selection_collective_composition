import numpy as np
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
def run(n):
    #parameter prepare
    mu=1e-4
    r=0.5
    s=3e-2
    N0=1000
    mbar=100
    ncomm=10
    rhat=0 if s>0 else 1
    nens=1000
    ncycle=10
    tcycle=np.log(n)/r
    #tcycle=np.log(100)/r
    
    fbins=np.linspace(0,1,100)
    
    #f0=0.1
    #frequency loop
    fsimdat=[]
    fsimrawdat=[]
    fselrawdat=[]
    fseldat=[]
    fextdat=[]
    fthdat=[]
    means=[]
    medians=[]
    mediansex=[]
    meansex=[]
    from tqdm import tqdm
    #for f0 in np.arange(0.00,0.99,0.01): 
    for f0 in tqdm(np.arange(0.01,0.99,0.01),leave=False):#[::int(np.sign(s))]): 
        fbins=np.linspace(np.maximum(0,f0-0.15),np.minimum(1,f0+0.15),30)
        #fbins=np.linspace(np.maximum(0,f0-0.05),np.minimum(1,f0+0.05),30)
        skewu=0.1 if s>0 else 0
        skewd=0.1 if s<0 else 0
        foribins=np.linspace(np.maximum(0,f0-0.15-skewd),np.minimum(1,f0+0.15+skewu),30)     #extend the range to describe the distribution more
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
        for e in range(nens):
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
    
        fsimrawdat.append(dats)
        simdat,_=np.histogram(dats,bins=foribins,density=True)
        #print('simulation histogram normalization test:',foribins,simdat)
        fsimdat.append(simdat)    #store histogram of simulation data
        cfs_sel=cf(w_sel,m_sel)
        fselrawdat.append(cfs_sel)    
        dat,_=np.histogram(cfs_sel,bins=fbins,density=True)
        fseldat.append(dat)    #store histogram of selected data
        medians.append(np.median(cfs_sel))
        means.append(np.mean(cfs_sel))
        ''' 
        ############################
        # Evaluation part
        ############################    
        pdfforplot=transition_probability_iid_pdf(foribins[:-1],f0,tcycle,N0,r,mu,s)
        #plt.plot(foribins[:-1],pdfforplot)
        #plt.hist(dats,bins=foribins,density=True)
        pdf=transition_probability_iid_pdf(fbins[:-1],f0,tcycle,N0,r,mu,s)
        pdfnorm=simpson(pdf,dx=fbins[1]-fbins[0])
        if f0-0.15<0 or 1<f0+0.15:
            pdf=pdf/pdfnorm        #if pdf range is not enough, reweight the pdf for sure completeness
        #pdf=transition_probability_iid_pdf(foribins[:-1],f0,tcycle,N0,r,mu,s)
        #pdfnorm=simpson(pdf,dx=fbins[1]-fbins[0])
        fthdat.append(pdfforplot)    
        cdf=np.cumsum(pdf*(fbins[1]-fbins[0]))
        extpdf=ncomm*(0.5+np.sign(s)*(0.5-cdf))**(ncomm-1)*pdf #Same dimension with dat
        #extpdf=ncomm*(cdf)**(ncomm-1)*pdf #Same dimension with dat
        #plt.plot(fbins[:-1],cdf)
        #plt.plot(fbins[:-1],extpdf)
        #print('Theory for',f0,cdf[-1],fbins,pdf,cdf,extpdf,np.sum(extpdf)*(fbins[1]-fbins[0]),np.sign(s),flush=True)
        #plt.show()
        #continue
        fextdat.append(extpdf)    
        median=quantile(0.5,fbins,extpdf,x0=np.min(fbins))
        mediansex.append(median)
        mean=np.sum(fbins[:-1]*extpdf*(fbins[1]-fbins[0]))
        meansex.append(mean)
        ''' 
        #exit()
    #folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
    #folder="data/cond/conditional_probability_taulog%s_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(n,N0,r,mu,s,ncomm,nens)
    folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix2"%(N0,r,mu,s,ncomm,nens)
    #foldere="data/cond/conditional_probability_taulog%s_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fix"%(n,ncomm,r,mu,s,ncomm,nens)
    print(folder)
    
    np.savetxt(folder+"_simraw.dat",fsimrawdat)
    np.savetxt(folder+"_sim.dat",fsimdat)
    np.savetxt(folder+"_sel.dat",fseldat)
    np.savetxt(folder+"_selraw.dat",fselrawdat)
    np.savetxt(folder+"_smedian.dat",medians)
    np.savetxt(folder+"_smean.dat",means)
    
    #np.savetxt(folder+"_th.dat",fthdat)
    #np.savetxt(folder+"_ext.dat",fextdat)
    #np.savetxt(folder+"_emedian.dat",mediansex)
    #np.savetxt(folder+"_emean.dat",meansex)

#for n in tqdm(range(4,21,1)):#[6,8,20,40,60,80,100]:
#    run(int(n))

run(11)
