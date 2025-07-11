import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.integrate as intg
from scipy.integrate import simpson
import scipy.special as spc
from analytic_results import *

###################
#Sampling distributions
##################
def sampling_distribution_pdf(f,sig,fo,psig,pmsig,g):
    #sample_mean=st.truncnorm.pdf(f,0,1,loc=fo,scale=psig/np.sqrt(g))
    sample_mean=st.norm.pdf(f,loc=fo,scale=pmsig)
    sample_dist=st.chi.pdf(sig,g-1,scale=psig/np.sqrt(g-1))

    return sample_mean*sample_dist
#Pg(ff|f,s)
def growth_distribution_pdf(f,sig,ff,t,N0,r,mu,s):
    m0=f*N0
    w0=N0-m0    
    msig0=sig*N0/np.sqrt(3+2*f**2)    
    wsig0=msig0
    #dist=st.truncnorm.pdf(ff,0,1,loc=barc_th(w0,m0,t,r,mu,s),scale=np.sqrt(sig2c_th(w0,wsig0**2,m0,msig0**2,t,r,mu,s)))
    dist=st.norm.pdf(ff,loc=barc_th(w0,m0,t,r,mu,s),scale=np.sqrt(sig2c_th(w0,wsig0**2,m0,msig0**2,t,r,mu,s)))
    return dist    

###################
# Distribution fmin of P(f>f_min)<1/g
##################

#nascent dirac delta function in probabilistic sence
def ndirac(x,eps=1e-2):
        return 1./eps if np.fabs(x)<eps/2 else 0

#Q(1/g) from given gaussian distribution (f0,sig0)
def freq_min_th(g,f0,sig0):
    return st.norm.ppf(1/g,loc=f0,scale=sig0)    

def integrand_g_inverse(f,sig,ff,t,N0,r,mu,s,pf,psig,pmsig,g):
    fmin=freq_min_th(g,f,sig)
    fbartau=barc_th_v3(fmin,t,r,mu,s)
    delta=ndirac(fbartau-ff)
    if delta==0:
        return 0
    else:
        sig_scale=sig/(psig+1e-16)
        exponent=0.5*((f-pf)/(pmsig+1e-16))**2+0.5*(sig_scale)**2*(g-1)
        #expo= 0 if exponent>30 else np.exp(-exponent)
        expo= np.exp(-exponent)
        coef=np.power(0.5*(g-1),0.5*(g-1))*np.power(sig_scale,g-2)/(np.pi*pmsig*psig*spc.gamma(0.5*(g-1))+1e-16)
    return coef*expo*delta
    
def transition_probability_g_inverse(ff,fsel,t,g,N0,Nf,r,mu,s):
    
    #from the selected value, generate reproduction distribution
    msel=barRm_th(fsel,N0=N0)
    wsel=N0-msel
    mselsig=sigRm_th(fsel,Nf,N0=N0)
    mselmeansig=sigRmbar_th(fsel,Nf,N0,g)
    wselsig=mselsig

    fpopmean=msel/N0
    fpopmeansig=np.sqrt(mselmeansig)/N0
    #fpopmean=barc_th(wsel,msel,0,r,mu,s)
    fpopstd=np.sqrt(mselsig)/N0
    #print("transition probability for ",ff," -> ",fpopmeansig,fpopstd/g)
    #fpopstd=np.sqrt(sig2c_th(wsel,wselsig,msel,mselsig,0,r,mu,s))
    #t=np.log(g+1)/r

    #val,err=intg.dblquad(integrand,0,1,lambda x:0,lambda x:1,args=(ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g))
    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
    df=fbins[1]-fbins[0]
    sbins=np.linspace(0,fpopstd*10*(g-1),100)
    ds=sbins[1]-sbins[0]
    fns=np.zeros((len(fbins),len(sbins)))
    for idx,f in enumerate(fbins):
        for idy,svar in enumerate(sbins):
            fns[idx,idy]=integrand_g_inverse(f,svar,ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g)
    
    val=np.sum(fns)*df*ds
    
    return val
def transition_probability_g_inverse_pdf(ff,fsel,t,g,N0,Nf,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_g_inverse(f,fsel,t,g,N0,Nf,r,mu,s))

    return np.array(rst)

###################
# Distribution of \int (Rep->Ext->Growth(determinisitic)) 
##################


#Extreme value of distribution Pg=g(1-F)^{g-1}P()
#P(fbar,sig) : sampling distribution
def integrand_ext_deterministic(f,ff,t,N0,r,mu,s,pf,psig,g):
    #growth - from initial sample mean and 'std' f,sig , get extreme value
    #Extrimal value distribution
    
    P=st.truncnorm.pdf(f,-pf/psig,(1-pf)/psig,loc=pf,scale=psig)    
    omF=1-st.truncnorm.cdf(f,-pf/psig,(1-pf)/psig,loc=pf,scale=psig)
    fg=barc_th_v2(f,N0,t,r,mu,s)

    delta=ndirac(fg-ff)
    if delta==0:
        return 0
    else:
        return delta*P*g*omF**(g-1)
def transition_probability_ext_deterministic(ff,fsel,t,g,N0,Nf,r,mu,s):
    
    #from the selected value, generate reproduction distribution
    msel=barRm_th(fsel,N0=N0)
    wsel=N0-msel
    mselsig=sigRm_th(fsel,Nf,N0=N0)
    mselmeansig=sigRmbar_th(fsel,Nf,N0,g)
    wselsig=mselsig

    fpopmean=msel/N0
    fpopmeansig=np.sqrt(mselmeansig)/N0
    #fpopmean=barc_th(wsel,msel,0,r,mu,s)
    fpopstd=np.sqrt(mselsig)/N0
    #print("transition probability for ",ff," -> ",fpopmeansig,fpopstd/g)
    #fpopstd=np.sqrt(sig2c_th(wsel,wselsig,msel,mselsig,0,r,mu,s))
    #t=np.log(g+1)/r


    #val,err=intg.dblquad(integrand,0,1,lambda x:0,lambda x:1,args=(ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g))
    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
    df=fbins[1]-fbins[0]
    #sbins=np.linspace(0,fpopstd*10*(g-1),100)
    #ds=sbins[1]-sbins[0]
    fns=np.zeros(len(fbins))
    for idx,f in enumerate(fbins):
        #for idy,svar in enumerate(sbins):
        fns[idx]=integrand_ext_deterministic(f,ff,t,N0,r,mu,s,fpopmean,fpopstd,g)
    
    val=np.sum(fns)*df
    
    return val
def transition_probability_ext_deterministic_pdf(ff,fsel,t,g,N0,Nf,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_ext_deterministic(f,fsel,t,g,N0,Nf,r,mu,s))

    return np.array(rst)

###################
# Distribution of \int (Rep->Ext->Growth(stochastic)) 
##################


#Extreme value of distribution Pg=g(1-F)^{g-1}P()
#P(fbar,sig) : sampling distribution
def integrand_ext_stochastic(f,ff,t,N0,r,mu,s,pf,psig,g):
    #growth - from initial sample mean and 'std' f,sig , get extreme value
    #Extrimal value distribution
    
    #Probability to choose f from N(Pf,psig)
    P=st.norm.pdf(f,loc=pf,scale=psig)    
    omF=1-st.norm.cdf(f,loc=pf,scale=psig)

    #probability to choose ff from N(fg,varg) after growth
    fg=barc_th_v2(f,N0,t,r,mu,s)
    varg=sig2c_th_v2(f,0,N0,t,r,mu,s)
    Pg=st.norm.pdf(ff,loc=fg,scale=np.sqrt(varg))

    return Pg*P*g*omF**(g-1)
def transition_probability_ext_stochastic(ff,fsel,t,g,N0,Nf,r,mu,s):
    
    #from the selected value, generate reproduction distribution
    msel=barRm_th(fsel,N0=N0)
    wsel=N0-msel
    mselsig=sigRm_th(fsel,Nf,N0=N0)
    mselmeansig=sigRmbar_th(fsel,Nf,N0,g)
    wselsig=mselsig

    fpopmean=msel/N0
    fpopmeansig=np.sqrt(mselmeansig)/N0
    #fpopmean=barc_th(wsel,msel,0,r,mu,s)
    fpopstd=np.sqrt(mselsig)/N0
    #print("transition probability for ",ff," -> ",fpopmeansig,fpopstd/g)
    #fpopstd=np.sqrt(sig2c_th(wsel,wselsig,msel,mselsig,0,r,mu,s))
    #t=np.log(g+1)/r


    #val,err=intg.dblquad(integrand,0,1,lambda x:0,lambda x:1,args=(ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g))
    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
    df=fbins[1]-fbins[0]
    #sbins=np.linspace(0,fpopstd*10*(g-1),100)
    #ds=sbins[1]-sbins[0]
    fns=np.zeros(len(fbins))
    for idx,f in enumerate(fbins):
        #for idy,svar in enumerate(sbins):
        fns[idx]=integrand_ext_stochastic(f,ff,t,N0,r,mu,s,fpopmean,fpopstd,g)
    
    val=np.sum(fns)*df
    
    return val
def transition_probability_ext_stochastic_pdf(ff,fsel,t,g,N0,Nf,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_ext_stochastic(f,fsel,t,g,N0,Nf,r,mu,s))

    return np.array(rst)

###################
# Distribution of \int (Rep(BN)->Ext->Growth(determinisitic)) 
##################


#Extreme value of distribution Pg=g(1-F)^{g-1}P()
#P(fbar,sig) : sampling distribution
def integrand_bin_ext_deterministic(f,ff,t,N0,r,mu,s,pf,psig,g):
    #growth - from initial sample mean and 'std' f,sig , get extreme value
    #Extrimal value distribution
    
    P=st.norm.pdf(f,loc=pf,scale=psig)    
    omF=1-st.norm.cdf(f,loc=pf,scale=psig)
    fg=barc_th_v2(f,N0,t,r,mu,s)

    delta=ndirac(fg-ff)
    if delta==0:
        return 0
    else:
        return delta*P*g*omF**(g-1)
def transition_probability_bin_ext_deterministic(ff,fsel,t,g,N0,Nf,r,mu,s):
    
    #from the selected value, generate reproduction distribution - binomial
    msel=N0*fsel
    wsel=N0-msel
    mselsig=N0*fsel*(1-fsel)
    mselmeansig=mselsig/np.sqrt(g)
    wselsig=mselsig

    fpopmean=msel/N0
    fpopmeansig=np.sqrt(mselmeansig)/N0
    #fpopmean=barc_th(wsel,msel,0,r,mu,s)
    fpopstd=np.sqrt(mselsig)/N0
    #print("transition probability for ",ff," -> ",fpopmeansig,fpopstd/g)
    #fpopstd=np.sqrt(sig2c_th(wsel,wselsig,msel,mselsig,0,r,mu,s))
    #t=np.log(g+1)/r


    #val,err=intg.dblquad(integrand,0,1,lambda x:0,lambda x:1,args=(ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g))
    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
    df=fbins[1]-fbins[0]
    #sbins=np.linspace(0,fpopstd*10*(g-1),100)
    #ds=sbins[1]-sbins[0]
    fns=np.zeros(len(fbins))
    for idx,f in enumerate(fbins):
        #for idy,svar in enumerate(sbins):
        fns[idx]=integrand_bin_ext_deterministic(f,ff,t,N0,r,mu,s,fpopmean,fpopstd,g)
    
    val=np.sum(fns)*df
    
    return val
def transition_probability_bin_ext_deterministic_pdf(ff,fsel,t,g,N0,Nf,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_bin_ext_deterministic(f,fsel,t,g,N0,Nf,r,mu,s))

    return np.array(rst)
###################
# Distribution of \int (Rep->Growth) 
##################

def integrand(f,sig,ff,t,N0,r,mu,s,pf,psig,pmsig,g):
    
    #growth - from initial sample mean and 'std' f,sig , grow in gaussian form
    m0=f*N0
    w0=N0-m0
    wsig2=(N0*sig)**2
    fg=barc_th(w0,m0,t,r,mu,s)
    sig2g=sig2c_th(w0,wsig2,m0,wsig2,-wsig2,t,r,mu,s)

    #Sampling distribution sampling distribution    
    sig_scale=sig/(psig+1e-16)
    
    #exponent=0.5*(ff-fg)**2/sig2g+0.5*g*g*((f-pf)/psig)**2+0.5*(sig_scale)**2*(g-1)
    exponent=0.5*(ff-fg)**2/(sig2g+1e-16)+0.5*((f-pf)/(pmsig+1e-16))**2+0.5*(sig_scale)**2*(g-1)
    #expo= 0 if exponent>30 else np.exp(-exponent)
    expo= np.exp(-exponent)
    #print(f,sig,exponent,expo)
    
    coef=np.power(0.5*(g-1),0.5*(g-1))*np.power(sig_scale,g-2)/(np.pi*np.sqrt(sig2g)*pmsig*psig*spc.gamma(0.5*(g-1))+1e-16)
    #coef=np.sqrt(g)*np.power(0.5*(g-1),0.5*(g-1))*np.power(sig_scale,g-2)/(np.pi*np.sqrt(sig2g)*psig**2*spc.gamma(0.5*(g-1)))

    #if(coef*expo>0):
    #print(f,sig,coef,exponent,expo)
    return coef*expo

def transition_probability(ff,fsel,t,g,N0,Nf,r,mu,s):
    
    #from the selected value, generate reproduction distribution
    msel=barRm_th(fsel,N0=N0)
    wsel=N0-msel
    mselsig=sigRm_th(fsel,Nf,N0=N0)
    mselmeansig=sigRmbar_th(fsel,Nf,N0,g)
    wselsig=mselsig

    fpopmean=msel/N0
    fpopmeansig=np.sqrt(mselmeansig)/N0
    #fpopmean=barc_th(wsel,msel,0,r,mu,s)
    fpopstd=np.sqrt(mselsig)/N0
    #print("transition probability for ",ff," -> ",fpopmeansig,fpopstd/g)
    #fpopstd=np.sqrt(sig2c_th(wsel,wselsig,msel,mselsig,0,r,mu,s))
    #t=np.log(g+1)/r
    #t=0

    #val,err=intg.dblquad(integrand,0,1,lambda x:0,lambda x:1,args=(ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g))
    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
    df=fbins[1]-fbins[0]
    sbins=np.linspace(0,fpopstd*10*(g-1),100)
    ds=sbins[1]-sbins[0]
    fns=np.zeros((len(fbins),len(sbins)))
    for idx,f in enumerate(fbins):
        for idy,svar in enumerate(sbins):
            fns[idx,idy]=integrand(f,svar,ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g)
    
    val=np.sum(fns)*df*ds
    
    return val

def transition_probability_pdf(ff,fsel,t,g,N0,Nf,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability(f,fsel,t,g,N0,Nf,r,mu,s))

    return np.array(rst)

###################
# Distribution of \int (Rep(BN)->Growth) with iid 
##################

def integrand_iid(f,ff,t,N0,r,mu,s,barf0,sig2f0):

    #Reproduction
    rep=st.norm.pdf(f,loc=barf0,scale=np.sqrt(sig2f0))
        
    #rep=st.norm.pdf(f,loc=barf0,scale=np.sqrt(barf0*(1-barf0)/N0))
        
    #growth - from initial sample mean and 'std' f,sig , grow in gaussian form
    fg=barc_th_v2(f,N0,t,r,mu,s)
    sig2g=sig2c_th_v2(f,0,N0,t,r,mu,s)
    grt=st.norm.pdf(ff,loc=fg,scale=np.sqrt(sig2g))
        
    return rep*grt

def integrand_iid_alt(f,ff,t,N0,r,mu,s,barf0,sig2f0):
    #for sharp rep fucntion
    #Reproduction
    rep=st.norm.pdf(f,loc=barf0,scale=np.sqrt(sig2f0))
    if np.isnan(rep):
        rep=np.ones(shape=len(f))/(f[1]-f[0])
        
    #rep=st.norm.pdf(f,loc=barf0,scale=np.sqrt(barf0*(1-barf0)/N0))
        
    #growth - from initial sample mean and 'std' f,sig , grow in gaussian form
    fg=barc_th_v2(f,N0,t,r,mu,s)
    sig2g=sig2c_th_v2(f,sig2f0,N0,t,r,mu,s)
    grt=st.norm.pdf(ff,loc=fg,scale=np.sqrt(sig2g))
        
    return rep*grt

def transition_probability_iid(ff,fsel,t,N0,r,mu,s):
    
    #from the selected value, generate reproduction distribution
    barf0=fsel
    sig2m0=N0*fsel*(1-fsel)
    #sig2f0=(1-2*fsel+2*fsel**2)*sig2m0/N0**2    
    sig2f0=sig2m0/N0**2    

    val,err=intg.quad(integrand_iid,0,1,args=(ff,t,N0,r,mu,s,barf0,sig2f0))
    if np.isnan(val)  :
        val,err=intg.quad(integrand_iid_alt,0,1,args=(ff,t,N0,r,mu,s,barf0,sig2f0+1e-16))
        '''
        #simple gubungujukbub
        #grid
        fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
        df=fbins[1]-fbins[0]
        sbins=np.linspace(0,fpopstd*10*(g-1),100)
        ds=sbins[1]-sbins[0]
        fns=np.zeros((len(fbins),len(sbins)))
        for idx,f in enumerate(fbins):
            fns[idx]=integrand_iid(fsel,ff,t,N0,r,mu,s,barf0,sig2f0)
        
        val=np.sum(fns)*df*ds
        '''
    return val

def transition_probability_iid_pdf(ff,fsel,t,N0,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_iid(f,fsel,t,N0,r,mu,s))

    return np.array(rst)

###################
# Distribution of \int (Rep->Growth->ext) 
##################
def integrand_ext(f,sig,ff,t,N0,r,mu,s,pf,psig,pmsig,g):
    
    #growth - from initial sample mean and 'std' f,sig , grow in gaussian form
    m0=f*N0
    w0=N0-m0
    wsig2=(N0*sig)**2
    fg=barc_th(w0,m0,t,r,mu,s)
    #sig2g=sig2c_th(w0,wsig2,m0,wsig2,-wsig2,t,r,mu,s)
    sig2g=sig2c_th_v2(f,sig,N0,t,r,mu,s)

    #Sampling distribution sampling distribution    
    sig_scale=sig/(psig+1e-16)
    
    exponent=0.5*(ff-fg)**2/(sig2g+1e-16)+0.5*((f-pf)/(pmsig+1e-16))**2+0.5*(sig_scale)**2*(g-1)
    #expo= 0 if exponent>30 else np.exp(-exponent)
    expo= np.exp(-exponent)
    #print(f,sig,exponent,expo)
    
    sfval=st.norm.sf(ff,loc=fg,scale=np.sqrt(sig2g+1e-16))
    coef=g*np.power(sfval,g-1)*np.power(0.5*(g-1),0.5*(g-1))*np.power(sig_scale,g-2)/(np.pi*np.sqrt(sig2g)*pmsig*psig*spc.gamma(0.5*(g-1))+1e-16)
    #coef=g*np.power(sfval,g-1)*np.sqrt(g)*np.power(0.5*(g-1),0.5*(g-1))*np.power(sig_scale,g-2)/(np.pi*np.sqrt(sig2g)*psig**2*spc.gamma(0.5*(g-1)))

    return coef*expo

def transition_probability_ext(ff,fsel,t,g,N0,Nf,r,mu,s):
    
    #from the selected value, generate reproduction distribution
    msel=barRm_th(fsel,N0=N0)
    wsel=N0-msel
    mselsig=sigRm_th(fsel,Nf,N0=N0)
    mselmeansig=sigRmbar_th(fsel,Nf,N0,g)
    wselsig=mselsig

    fpopmean=msel/N0
    fpopmeansig=np.sqrt(mselmeansig)/N0
    #print("transition probability ext for ",ff," -> ",mselmeansig)
    #fpopmean=barc_th(wsel,msel,0,r,mu,s)
    fpopstd=np.sqrt(mselsig)/N0
    #fpopstd=np.sqrt(sig2c_th(wsel,wselsig,msel,mselsig,0,r,mu,s))
    #t=np.log(g+1)/r
    
    #val,err=intg.dblquad(integrand_ext,0,1,lambda x:0,lambda x:1,args=(ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g))
    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100)
    df=fbins[1]-fbins[0]
    sbins=np.linspace(0,fpopstd*10*(g-1),100)
    ds=sbins[1]-sbins[0]
    fns=np.zeros((len(fbins),len(sbins)))
    for idx,f in enumerate(fbins):
        for idy,svar in enumerate(sbins):
            fns[idx,idy]=integrand_ext(f,svar,ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g)
    
    val=np.sum(fns)*df*ds

    return val
def transition_probability_ext_pdf(ff,fsel,t,g,N0,Nf,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_ext(f,fsel,t,g,N0,Nf,r,mu,s))

###################
# Distribution of pdf \int (Rep->Growth->ext) : Gaussian approximation, Not extreme value

#integrand : 1/(1-f)^2 * (wbar(t)*sig2m(t)+f/(1-f)mbar(t)*sig2w(t))/(sig2m+(f/(1-f))^2sig2w(t)) * Nm'(m|mb'(t),sigm'(t))

#integrator : dm, dm0

##################
def integrand_gaussian(m0,ff,t,N0,r,mu,s,barm0,sig2m0):

    barw=barw_th(N0-m0,t,r,mu)
    sig2w=sig2w_th(N0-m0,0,t,r,mu)
    barm=barm_th(m0,N0-m0,t,r,mu,s)
    sig2m=sig2m_th(m0,0,N0-m0,0,0,t,r,mu,s)

    fp1mf=ff/(1-ff+1e-16) #prevent diverge
    
    coef=(barw*sig2m+fp1mf*barm*sig2w)/((1-ff+1e-16)**2*(sig2m+fp1mf**2*sig2w))    
    normal=st.norm.pdf(m0,loc=barm0,scale=np.sqrt(sig2m0))
    
    return coef*normal

def transition_probability_gaussian(ff,fsel,t,N0,r,mu,s):

    #For dm0
    barm0=N0*fsel
    sig2m0=N0*fsel*(1-fsel)
    print(N0,fsel,sig2m0)
    
    #Gaussian quardrature    
    val,err=intg.quad(integrand_gaussian,np.maximum(barm0-3*sig2m0,0),barm0+3*sig2m0,args=(ff,t,N0,r,mu,s,barm0,sig2m0))
    '''    
    #simple gubungujukbub
    #grid
    fbins=np.linspace(np.maximum(0,fsel-5*fpopstd),np.minimum(1,fsel+5*fpopstd),100
    df=fbins[1]-fbins[0]
    sbins=np.linspace(0,fpopstd*10*(g-1),100)
    ds=sbins[1]-sbins[0]
    fns=np.zeros((len(fbins),len(sbins)))
    for idx,f in enumerate(fbins):
        for idy,svar in enumerate(sbins):
            fns[idx,idy]=integrand_ext(f,svar,ff,t,N0,r,mu,s,fpopmean,fpopstd,fpopmeansig,g)
    
    val=np.sum(fns)*df*ds
    '''

    return val
def transition_probability_gaussian_pdf(ff,fsel,t,N0,r,mu,s):
    rst=[]
    for f in ff:
        rst.append(transition_probability_gaussian(f,fsel,t,N0,r,mu,s))

    return np.array(rst)


####
# ETC
####
def transition_probability_cdf_from_data(ff,fsel,t,g,N0,Nf,r,mu,s):

    pdfs=transition_probability_cdf_from_data(ff,fsel,t,g,N0,Nf,r,mu,s)
    cdfs=np.cumsum(pdfs)*(ff[1:]-ff[:-1])
    return cdfs,pdfs

def integrand2(ff,fhat,fsig,fsel,t,g,N0,Nf,r,mu,s):
    score=st.norm.pdf(ff,loc=fhat,scale=fsig)
    tr=transition_probability(ff,fsel,t,g,N0,Nf,r,mu,s)
    return score*tr

def transition_probability_extreme_min(cdf,pdf,g):
    return g*np.power(np.fabs(1-cdf),g-1)*pdf

#print("checkit",intg.quad(integrand2,0,1,args=(0.1,0.1,0.1,ncomm,N0,(ncomm+1)*N0,r,mu,s)))

#sampler with gaussian function
def selection_weight(f,fhat,fsig):    
    return st.norm.pdf(f,loc=fhat,scale=fsig)


if __name__=='__main__':
    #parameter prepare
    mu=1e-4
    r=0.5
    s=-3e-2
    N0=1000
    mbar=100
    ncomm=10
    g=ncomm
    rhat=0
    ncycle=10
    t=np.log(ncomm+1)/r

    def freq(w,m):
        return np.divide(m,w+m)
    
    ############################
    # Simulation part
    ############################
    from model import model
    from tau_leaping import tau_leaping
    from selection import select
    from selection import community_function as cf
    from selection import score_function as sf
    from growth_alter import *
    
    #model
    proFunc=[
    lambda x: r*x[0],
    lambda x: mu*x[0],
    lambda x: (r+s)*x[1]
    ]
    changeVec=np.array([[1,0],[-1,1],[0,1]]).astype(int)
    rxnOrder=np.array([[1,0],[1,0],[0,1]]).astype(int)
    growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)
    solver=tau_leaping(model=growthmodel,eps=0.03)
    
    #Start from give initial state
    fsel=0.20
    nens=1000
    #For save all data
    lastw=np.array([])
    lastm=np.array([])
    minf=[]
    from tqdm import tqdm
    for e in tqdm(range(nens)):
    #for e in range(nens):
        #print("ensemble ",e)
        
        #reproduction phase    
        m0sel=np.random.binomial(N0,fsel,size=ncomm)
        w0sel=N0-m0sel
    
        #growth phase
        for j in range(ncomm):    
            #T,X=solver.run(np.array([w0sel[j],m0sel[j]]),0,t)
            wsamp=growth_sampling_w(w0sel[j],0,t,r,mu,s)
            msamp=growth_sampling_m(m0sel[j],0,w0sel[j],0,0,t,r,mu,s)
            #print(X[0,-1]+X[1,-1])
            #wsel=np.append(wsel,X[0,:-1])
            #msel=np.append(msel,X[1,:-1])
        
            #store for making distribution
            lastw=np.append(lastw,wsamp)#X[0,-1])
            lastm=np.append(lastm,msamp)#X[1,-1])    

        ind_sel=select(lastw[-ncomm:],lastm[-ncomm:],r=r,s=s,rhat=0)
        minf.append(cf(lastw[ind_sel+e*ncomm],lastm[ind_sel+e*ncomm]))
    cfs=cf(lastw,lastm) #get frequency
    cfmean=np.mean(cfs) #get mean
    cfstd=np.std(cfs)    #get std    

    fbins=np.linspace(np.min(cfs)-0.05,np.max(cfs)+0.05,30)
    fbinsex=np.linspace(np.min(minf)-0.05,np.max(minf)+0.05,30)
    fbinmid=fbins[:-1]+0.5*(fbins[1]-fbins[0])
    fbinexmid=fbinsex[:-1]+0.5*(fbinsex[1]-fbinsex[0])
    fbinsfine=np.linspace(np.min(cfs)-0.05,np.max(cfs)+0.05,100)
    fbinsfinemid=fbinsfine[:-1]+0.5*(fbinsfine[1]-fbinsfine[0])
    dat,_,_=plt.hist(cfs,bins=fbins,density=True,alpha=0.5,histtype='step',label='Simulation')
    plt.plot(fbinmid,st.norm.pdf(fbinmid,loc=np.mean(cfs),scale=np.std(cfs)),color='green',lw=1,label='Fitting') # fitting from stochastic data
    
    pdf=transition_probability_iid_pdf(fbinmid,fsel,t,N0,r,mu,s)
    plt.plot(fbinmid,pdf,c='black',lw=1,label='Theory')
    
    barc=barc_th_v2(fsel,N0,t,r,mu,s)
    sig2c=sig2c_th_v2(fsel,fsel*(1-fsel)/N0,N0,t,r,mu,s)
    print('mean=',np.mean(cfs),barc,'std=',np.std(cfs),np.sqrt(sig2c))
    
    datex,_,_=plt.hist(minf,bins=fbinsex,density=True,alpha=0.5,histtype='step',label='Simex')
    
    pdf=transition_probability_iid_pdf(fbinsfinemid,fsel,t,N0,r,mu,s)
    pdfnorm=simpson(pdf,dx=fbinsfine[1]-fbinsfine[0])
    print(pdfnorm)
    pdf=pdf/pdfnorm
    cdf=np.cumsum(pdf*(fbinsfine[1]-fbinsfine[0]))
    print(cdf[-1])
    cdf=cdf/cdf[-1]
    extpdf=ncomm*(1-cdf)**(ncomm-1)*pdf #Same dimension with dat
    plt.plot(fbinsfinemid,extpdf,c='red',lw=1,label='Theo_ext')

    from stats_custom_tools import * 
    extmedian=quantile(0.5,fbinsfine,extpdf,x0=np.min(fbinsfine))
    extmean=simpson(fbinsfine[:-1]*extpdf,dx=fbinsfine[1]-fbinsfine[0])

    print(np.mean(minf),extmean,np.median(minf),extmedian)    
    plt.legend()
    
    plt.show()
